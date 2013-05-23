from __future__ import print_function

# symodesys imports
from symodesys.integrator import IVP_Integrator
from symodesys.helpers import import_

# stdlib imports
import tempfile
import shutil
import re
import os
from collections import OrderedDict
from functools import reduce
from operator import add

# other imports
import sympy
import numpy as np
import mako
import mako.template

from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import Extension
from Cython.Distutils import build_ext
import cython_gsl

class Generic_Code(object):
    """
    Wraps some sympy functionality of code generation from matrices
    returned by FirstOrderODESystem.dydt and FirstOrderODESystem.dydt_jac
    """

    def __init__(self, fo_odesys, tempdir = None, save_temp = False):
        self._fo_odesys = fo_odesys
        if tempdir:
            self._tempdir = tempdir
        else:
            self._tempdir = tempfile.mkdtemp()
            self._remove_tempdir_on_clean = True
        self._save_temp = save_temp

        if not os.path.isdir(self._tempdir):
            os.makedirs(self._tempdir)
            self._remove_tempdir_on_clean = True
        self._written_files = []
        self._source_files = []
        self._write_code()

    def _write_code(self):
        for path in self.copy_files:
            srcpath = os.path.join(os.path.dirname(__file__), path)
            dstpath = os.path.join(self._tempdir,
                         os.path.basename(path).replace('_template', ''))
            shutil.copy(srcpath, dstpath)
            if path in self._ori_sources:
                self._source_files.append(dstpath)

        for k, (path, attr) in self.templates.iteritems():
            srcpath = os.path.join(os.path.dirname(__file__), path)
            outpath = os.path.join(self._tempdir,
                         os.path.basename(path).replace('_template', ''))
            subs = getattr(self, attr)
            template = mako.template.Template(open(srcpath, 'rt').read())
            try:
                rendered = template.render(**subs)
            except:
                print(mako.exceptions.text_error_template().render())
                raise
            open(outpath, 'wt').write(rendered)
            if path in self._ori_sources:
                self._source_files.append(outpath)
            self._written_files.append(outpath)

    def compile_and_import_binary(self):
        binary_path = self._compile()
        return import_(binary_path)

    def _compile(self):
        extension_name = 'pyinterface'
        setup(
            script_name =  'DUMMY_SCRIPT_NAME',
            script_args =  ['build_ext',  '--build-lib', self._tempdir],
            include_dirs = [cython_gsl.get_include()],
            cmdclass = {'build_ext': build_ext},
            ext_modules = [
                Extension(
                    extension_name,
                    self._source_files,
                    libraries=cython_gsl.get_libraries(),
                    library_dirs=[cython_gsl.get_library_dir()],
                    include_dirs=[cython_gsl.get_cython_include_dir()]),
                ]
            )
        return os.path.join(self._tempdir, extension_name)

    def clean(self):
        if not self._save_temp:
            map(os.unlink, self._written_files)
            if self._remove_tempdir_on_clean:
                shutil.rmtree(self._tempdir)

    def __del__(self):
        """
        When Generic_Code object is collected by GC
        self._tempdir is deleted
        """
        self.clean()

    @property
    def NY(self):
        return len(self._fo_odesys.non_analytic_depv)

    @property
    def cprog_param_symbs(self):
        """
        Bridge slight difference in meaning of parameter going from
        ODESystem to numerical integration routine (analytically
        solving part of ODESys introduces new 'parameters')
        """
        return self._fo_odesys.param_and_sol_symbs

    @property
    def ccode_func(self):
        non_analytic_expr = self._fo_odesys.non_analytic_f.values()
        cse_defs, cse_exprs = sympy.cse(
            non_analytic_expr, symbols = sympy.numbered_symbols('cse'))

        f_cexprs = [self.arrayify(sympy.ccode(x)) for x in cse_exprs]

        cse_func = []
        for var_name, var_expr in cse_defs:
            c_var_expr = self.arrayify(sympy.ccode(var_expr))
            cse_func.append((var_name, c_var_expr))

        return {'cse_func': cse_func, 'f': f_cexprs}

    @property
    def ccode_ode(self):
        return {'NY': self.NY}

    @property
    def ccode_main_ex(self):
        params = [(str(p), 1.0) for p in self.cprog_param_symbs]
        y0 = ', '.join(['1.0'] * self.NY)
        params = ', '.join(['1.0'] * len(self.cprog_param_symbs))
        return {'NY': self.NY,
                'Y0_COMMA_SEP_STR': y0,
                'PARAM_VALS_COMMA_SEP_STR': params}

    @property
    def ccode_jac(self):
        na_jac = self._fo_odesys.non_analytic_jac
        na_f = map(self._fo_odesys.unfunc_depv, self._fo_odesys.non_analytic_f.values())
        indepv = self._fo_odesys.indepv

        sparse_jac = OrderedDict(reduce(add, [
            [((i, j), expr) for j, expr in enumerate(row) if expr != 0]\
            for i, row in enumerate(na_jac.tolist())
            ]))
        dfdt = OrderedDict([
            (i, expr) for i, expr in enumerate(
                [x.diff(indepv) for x in na_f]
                ) if expr != 0
            ])
        cse_defs, cse_exprs = sympy.cse(
            sparse_jac.values() + dfdt.values(),
            symbols = sympy.numbered_symbols('cse')
            )

        jac_cexprs = zip(sparse_jac.keys(), [
            self.arrayify(sympy.ccode(x)) for x \
            in cse_exprs[:len(sparse_jac)]
            ])

        dfdt_cexprs = zip(dfdt.keys(), [
            self.arrayify(sympy.ccode(x)) for x \
            in cse_exprs[len(sparse_jac):]
            ])

        cse_jac = []
        for var_name, var_expr in cse_defs:
            c_var_expr = self.arrayify(sympy.ccode(var_expr))
            cse_jac.append((var_name, c_var_expr))

        return {'cse_jac': cse_jac, 'jac': jac_cexprs,
                'dfdt': dfdt_cexprs,
                'NY': self.NY}


    def arrayify(self, ccode):
        for i, depv in enumerate(self._fo_odesys.non_analytic_depv):
            ccode = ccode.replace(str(depv), 'y[{}]'.format(i))
        for i, depv in enumerate(self._fo_odesys.analytic_depv):
            ccode = ccode.replace(
                str(depv), depv.func.__name__ + '(t, y, params)')
        for i, cparam in enumerate(self.cprog_param_symbs):
            ccode = ccode.replace(str(cparam), 'k[{}]'.format(i))
        ccode = re.sub('_(\d+)', r'[\1]', ccode)
        return ccode


class GSL_Code(Generic_Code):

    # Implement hash of fo_odesys and hash of code?
    # Serialization to double check against collision?


    copy_files = ('drivers.c', 'pyinterface.pyx') + \
                 ('drivers.h', 'func.h', 'jac.h', 'Makefile')

    templates = {
        'dydt': ('func_template.c', 'ccode_func'),
        'dydt_jac': ('jac_template.c', 'ccode_jac'),
        'main_ex': ('main_ex_template.c', 'ccode_main_ex')
        }

    _ori_sources = list(copy_files[:2]) + [templates[k][0] for k in ['dydt', 'dydt_jac']]



class GSL_IVP_Integrator(IVP_Integrator):
    """
    IVP integrator using GNU Scientific Library routines odeiv2

    remember to run on the instance on program exit.
    _code.clean()
    """

    # step_type choices are in `step_types` in pyinterface.pyx
    integrate_args = {'step_type': (
        'rk2','rk4','rkf45','rkck','rk8pd','rk2imp',
        'rk4imp','bsimp','rk1imp','msadams','msbdf'),
    }

    def __init__(self, **kwargs):
        self.tempdir = kwargs.pop('tempdir', None)
        self.save_temp = kwargs.pop('save_temp', False)
        super(GSL_IVP_Integrator, self).__init__(**kwargs)

    def set_fo_odesys(self, fo_odesys):
        super(GSL_IVP_Integrator, self).set_fo_odesys(fo_odesys)
        self._binary = None # <-- Clears cache
        self._code = GSL_Code(self._fo_odesys,
                              tempdir = self.tempdir,
                              save_temp = self.save_temp)

    @property
    def binary(self):
        if self._binary == None:
            self._binary = self._code.compile_and_import_binary()
        return self._binary


    def run(self, y0, t0, tend, param_vals, N,
                  abstol = None, reltol = None, h = None, h_init = None, h_max=0.0, order = 0):
        """
        hmax won't be set if 0.0
        """
        if h_init == None:
            h_init = 1e-9 # TODO: along the lines of: h_init=f(y0, dydt, jac, abstol, reltol)

        # Set (c-program) init vals to the possible subset of non_anlytic depv
        y0_arr = np.array(
            [y0[k] for k in self._fo_odesys.non_analytic_depv],
            dtype = np.float64)

        # Extend (c-program) params with y0 values of analytic functions
        cprog_param_vals = dict(
            param_vals.items() + \
            [(k, v) for k, v in y0.items() if k \
             in self._fo_odesys.analytic_depv])

        # C-program knows nothing about dicts, provide values as an array
        params_arr = np.array([cprog_param_vals[k] for k \
                               in self._code.cprog_param_symbs],
                              dtype = np.float64)
        if N > 0:
            # Fixed stepsize
            self.init_Yout_tout_for_fixed_step_size(t0, tend, N, order)
            tout, Yout = self.binary.integrate_equidistant_output(
                t0, tend, y0_arr, N, h_init, h_max, self.abstol, self.reltol, params_arr, order)
            self.tout = tout
            self.Yout = Yout
        else:
            raise NotImplementedError
