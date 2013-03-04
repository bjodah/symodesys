
# symodesys imports
from symodesys.integrator import IVP_Integrator
from symodesys.helpers import import_

# stdlib imports
import tempfile
from shutil import rmtree

# other imports
from mako.template import Template


from distutils.core import setup
from Cython.Distutils import Extension
from Cython.Distutils import build_ext
import cython_gsl

class Genric_Code(object):
    """
    Wraps some sympy functionality of code generation from matrices
    returned by FirstOrderODESystem.dydt and FirstOrderODESystem.dydt_jac
    """

    # CSE


class GSL_Code(object):

    # Implement hash of fo_odesys and hash of code?
    # Serialization to double check against collision?

    templates = {'dydt': 'gsl/func_template.c',
                 'dydt_jac': 'gsl/jac_template.c',
                 'ode': 'ode_template.c'}

    subs = {'dydt': ['f', 'cses'],
            'dydt_jac': ['jac', 'dfdt', 'NY'],
            'ode': ['NY']}

    @property
    def NY(self):
        return len(self.fo_odesys.dep_var_func_symbs)

    @property
    def dydt(self):
        pass

    @property
    def cse_func(self):
        pass

    @property
    def jac(self):
        pass

    @property
    def dfdt(self):
        pass

    @property
    def cse_jac(self):
        pass

    def __init__(self, fo_odesys, tempdir = None):
        self._fo_odesys = fo_odesys
        self._tempdir = tempdir or tempfile.mkdtemp()
        assert os.path.isdir(self._tempdir)
        self._write_code()

    def _write_code(self):
        for k, v in self.templates.iteritems():
            srcpath = os.path.join(os.path.dirname(self.__file__), v)
            outpath = os.path.join(self._tempdir, os.paht.basename(v))
            subs = {s: getattr(self, s) for s in self.subs[k]}
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
            ext_modules = [Extension("odeiv",
                                     ["odeiv.pyx"],
                                     libraries=cython_gsl.get_libraries(),
                                     library_dirs=[cython_gsl.get_library_dir()],
                                     include_dirs=[cython_gsl.get_cython_include_dir()]),
                           Extension("OdeSystem",
                                     ["OdeSystem.pyx"],
                                     libraries=cython_gsl.get_libraries(),
                                     library_dirs=[cython_gsl.get_library_dir()],
                                     include_dirs=[cython_gsl.get_cython_include_dir()])]
            )

    def clean(self):
        rmtree(self._tempdir)

    def __del__(self):
        """ When Generic_Code object is collected by GC self._tempdir is deleted """
        self.clean()


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
