
# symodesys imports
from symodesys.integrator import IVP_Integrator
from symodesys.helpers import import_

# stdlib imports
import tempfile
from shutil import rmtree
import re
import os

# other imports
import sympy
from mako.template import Template


from distutils.core import setup
from Cython.Distutils import Extension
from Cython.Distutils import build_ext
import cython_gsl

class Generic_Code(object):
    """
    Wraps some sympy functionality of code generation from matrices
    returned by FirstOrderODESystem.dydt and FirstOrderODESystem.dydt_jac
    """

    def __init__(self, fo_odesys, tempdir = None):
        self._fo_odesys = fo_odesys
        self._tempdir = tempdir or tempfile.mkdtemp()
        assert os.path.isdir(self._tempdir)
        self._write_code()

    def _write_code(self):
        for k, (path, attr) in self.templates.iteritems():
            srcpath = os.path.join(os.path.dirname(__file__), path)
            outpath = os.path.join(self._tempdir, os.path.basename(path))
            subs = getattr(self, attr)
            template = Template(open(srcpath, 'rt').read())
            open(outpath, 'wt').write(template.render(**subs))

    def compile_and_import_binary(self):
        binary_path = self._compile()
        return import_(binary_path)

    def _compile(self):
        setup(
            script_name =  'DUMMY_SCRIPT_NAME',
            script_args =  ['build_ext',  '--inplace'],
            include_dirs = [cython_gsl.get_include()],
            cmdclass = {'build_ext': build_ext},
            ext_modules = [
                Extension(
                    "odeiv",
                    ["odeiv.pyx"],
                    libraries=cython_gsl.get_libraries(),
                    library_dirs=[cython_gsl.get_library_dir()],
                    include_dirs=[cython_gsl.get_cython_include_dir()]),
                Extension(
                    "OdeSystem",
                    ["OdeSystem.pyx"],
                    libraries=cython_gsl.get_libraries(),
                    library_dirs=[cython_gsl.get_library_dir()],
                    include_dirs=[cython_gsl.get_cython_include_dir()])
                ]
            )

    def clean(self):
        rmtree(self._tempdir)

    def __del__(self):
        """
        When Generic_Code object is collected by GC
        self._tempdir is deleted
        """
        self.clean()

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
            non_analytic_expr, symbols = sympy.numbered_symbols('cse_'))

        f_cexprs = [self.arrayify(sympy.ccode(x)) for x in cse_exprs]

        cse_func = []
        for var_name, var_expr in cse_defs:
            c_var_expr = self.arrayify(sympy.ccode(var_expr))
            cse_func.append((var_name, c_var_expr))

        return {'cse_func': cse_func, 'f': f_cexprs}

    @property
    def ccode_ode(self):
        return {'NY': len(self._fo_odesys.non_analytic_depv)}

    @property
    def ccode_jac(self):
        na_jac = self._fo_odesys.non_analytic_jacobian
        na_f = self._fo_odesys.non_analytic_f

        sparse_jac = OrderedDict([
            ((i, j), expr) for j, expr in enumerate(row) for i, row in \
            enumerate(na_jac.tolist()) if expr != 0])
        dfdt = OrderedDict([
            (i, expr) for i, expr in enumerate(
                [x.diff(indepv) for x in na_f]
                ) if expr != 0
            ])
        cse_defs, cse_exprs = sympy.cse(
            sparse_jac.values() + dfdt.values(),
            symbols = sympy.numbered_symbols('cse_')
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
                'NY': len(self._fo_odesys.non_analytic_depv)}


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

    templates = {'dydt': ('gsl/func_template.c', 'ccode_func'),
                 'dydt_jac': ('gsl/jac_template.c', 'ccode_jac'),
                 'ode': ('ode_template.c', 'ccode_ode')
                 }


class GSL_IVP_Integrator(IVP_Integrator):
    """
    IVP integrator using GNU Scientific Library routines odeiv2
    """

    def post_init(self):
        self._code = GSL_Code(self._fo_odesys)
        self._binary = self._code.compile_and_import_binary()

    def integrate(self, y0, t0, tend, N, h = None, order = 0):
        y0_val_lst = [y0[k] for k in self._fo_odesys.dep_var_func_symbs]
        if N > 0:
            # Fixed stepsize
            self.init_yout_tout_for_fixed_step_size(t0, tend, N, order)
            # Order give (super)dimensionality of yout
            h_init = 1e-10 # TODO: find h: max(dydt) = abstol
            h_max = 0.0 # hmax won't be set if 0.0
            print_values = False
            Yout = self._binary.integrate_ode_using_driver_fixed_step(
                t0, tend, y0_val_lst, N, param_lst, self.abstol, self.reltol,
                h_init, h_max, print_values)
        else:
            raise NotImplementedError
