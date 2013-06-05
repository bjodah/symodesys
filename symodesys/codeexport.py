from __future__ import print_function, division

# stdlib imports
import tempfile
import shutil
import re
import os
from collections import OrderedDict
from functools import reduce
from operator import add

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# External imports
import sympy
import mako
import mako.template

# Intrapackage imports
from symodesys.helpers import import_
from symodesys.integrator import IVP_Integrator


class Generic_Code(object):
    """
    Wraps some sympy functionality of code generation from matrices
    returned by FirstOrderODESystem.dydt and FirstOrderODESystem.dydt_jac
    """

    wcode = staticmethod(sympy.ccode)
    syntax = 'C'

    _basedir = '.' # source file paths relative to _basedir

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
        self._include_dirs = []
        self._libraries = []
        self._library_dirs = []
        self._include_dirs = []
        self._write_code()

    def _write_code(self):
        for path in self.copy_files:
            srcpath = os.path.join(self._basedir, path)
            dstpath = os.path.join(self._tempdir,
                         os.path.basename(path).replace('_template', ''))
            shutil.copy(srcpath, dstpath)
            if path in self._ori_sources:
                self._source_files.append(dstpath)

        for k, (path, attr) in self.templates.iteritems():
            srcpath = os.path.join(self._basedir, path)
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

    def _compile(self, extension_name = 'pyinterface'):
        setup(
            script_name =  'DUMMY_SCRIPT_NAME',
            script_args =  ['build_ext',  '--build-lib', self._tempdir],
            include_dirs = self._include_dirs,
            cmdclass = {'build_ext': build_ext},
            ext_modules = [
                Extension(
                    extension_name,
                    self._source_files,
                    libraries=self._libraries,
                    library_dirs=self._library_dirs,
                    include_dirs=self._include_dirs),
                ]
            )
        return os.path.join(self._tempdir, extension_name)

    def clean(self):
        """ Delete temp dir if not save_temp set at __init__ """
        if not self._save_temp:
            map(os.unlink, self._written_files)
            if self._remove_tempdir_on_clean:
                shutil.rmtree(self._tempdir)

    def __del__(self):
        """
        When Generic_Code object is collected by GC
        self._tempdir is (possibly) deleted
        """
        self.clean()

    @property
    def NY(self):
        return len(self._fo_odesys.non_analytic_depv)

    @property
    def prog_param_symbs(self):
        """
        Bridge slight difference in meaning of parameter going from
        ODESystem to numerical integration routine (analytically
        solving part of ODESys introduces new 'parameters')
        """
        return self._fo_odesys.param_and_sol_symbs

    @property
    def code_func(self):
        non_analytic_expr = self._fo_odesys.non_analytic_f.values()
        cse_defs, cse_exprs = sympy.cse(
            non_analytic_expr, symbols = sympy.numbered_symbols('cse'))

        f_cexprs = [self.arrayify(self.wcode(x)) for x in cse_exprs]

        cse_func = []
        for var_name, var_expr in cse_defs:
            c_var_expr = self.arrayify(self.wcode(var_expr))
            cse_func.append((var_name, c_var_expr))

        return {'cse_func': cse_func, 'f': f_cexprs}

    @property
    def code_ode(self):
        return {'NY': self.NY}

    @property
    def code_main_ex(self):
        params = [(str(p), 1.0) for p in self.prog_param_symbs]
        y0 = ', '.join(['1.0'] * self.NY)
        params = ', '.join(['1.0'] * len(self.prog_param_symbs))
        return {'NY': self.NY,
                'Y0_COMMA_SEP_STR': y0,
                'PARAM_VALS_COMMA_SEP_STR': params}

    @property
    def code_jac(self):
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
            self.arrayify(self.wcode(x)) for x \
            in cse_exprs[:len(sparse_jac)]
            ])

        dfdt_cexprs = zip(dfdt.keys(), [
            self.arrayify(self.wcode(x)) for x \
            in cse_exprs[len(sparse_jac):]
            ])

        cse_jac = []
        for var_name, var_expr in cse_defs:
            c_var_expr = self.arrayify(self.wcode(var_expr))
            cse_jac.append((var_name, c_var_expr))

        return {'cse_jac': cse_jac, 'jac': jac_cexprs,
                'dfdt': dfdt_cexprs,
                'NY': self.NY}


    def arrayify(self, scode):
        """
        Returns arrayified expression
        self.syntax='C' implies C syntax, 'F' fortran respectively.
        """
        for i, depv in enumerate(self._fo_odesys.non_analytic_depv):
            tgt = {'C':'y[{}]', 'F':'y({})'}.get(self.syntax)
            scode = scode.replace(str(depv), tgt.format(i))
        # for i, depv in enumerate(self._fo_odesys.analytic_depv):
        #     scode = scode.replace(
        #         str(depv), depv.func.__name__ + '(t, y, params)')
        for i, param in enumerate(self.prog_param_symbs):
            tgt = {'C':'k[{}]', 'F':'k({})'}.get(self.syntax)
            scode = scode.replace(str(param), tgt.format(i))
        tgt = {'C':r'[\1]', 'F':r'(\1)'}.get(self.syntax)
        scode = re.sub('_(\d+)', tgt, scode)
        return scode

class Binary_IVP_Integrator(IVP_Integrator):
    """
    Generic baseclass for writing integrators
    which generate code (set Generic_Code derived
    class as CodeClass in derived class).
    """

    CodeClass = None

    def __init__(self, **kwargs):
        self.tempdir = kwargs.pop('tempdir', None)
        self.save_temp = kwargs.pop('save_temp', False)
        super(Binary_IVP_Integrator, self).__init__(**kwargs)

    def set_fo_odesys(self, fo_odesys):
        super(Binary_IVP_Integrator, self).set_fo_odesys(fo_odesys)
        self._binary = None # <-- Clears cache
        self._code = self.CodeClass(self._fo_odesys,
                              tempdir = self.tempdir,
                              save_temp = self.save_temp)

    def clean(self):
        self._code.clean()

    @property
    def binary(self):
        """
        Returns compiled and imported module.
        Note: lazy caching is employed, set self._binary equal to None to invalidate
        """
        if self._binary == None:
            self._binary = self._code.compile_and_import_binary()
        return self._binary
