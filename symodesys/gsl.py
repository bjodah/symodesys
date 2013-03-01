
# symodesys imports
from symodesys.integrator import IVP_Integrator

# stdlib imports
import tempfile
from shutil import rmtree

# other imports
import mako

class Genric_Code(object):
    """
    Wraps some sympy functionality of code generation from matrices
    returned by FirstOrderODESystem.dydt and FirstOrderODESystem.dydt_jac
    """

    # CSE


class GSL_Code(object):

    # Implement hash of fo_odesys and hash of code?
    # Serialization to double check against collision?

    templates = {'dydt' = 'gsl/func_template.c',
                 'dydt_jac': 'gsl/jac_template.c',
                 'ode': 'ode_template.c'}

    subs = {'dydt' = ['f', 'cses'],
            'dydt_jac': ['jac', 'dfdt', 'NY'],
            'ode': ['NY']}

    @Property
    def NY(self):
        return len(self.fo_odesys.dep_var_func_symbs)

    @Property
    def dydt(self):
        pass

    @Property
    def cse_func(self):
        pass

    @Property
    def jac(self):
        pass

    @Property
    def dfdt(self):
        pass

    @Property
    def cse_jac(self):
        pass


    def __init__(self, fo_odesys, tempdir = None):
        self._fo_odesys = fo_odesys
        self._tempdir = tempdir or tempfile.mkdtemp()
        assert os.path.isdir(self._tempdir)
        self._generate_code()
        self._write_code()

    def clean(self):
        rmtree(self._tempdir)

    def compile_and_import_binary(self):
        self._compile()

    def _generate_code()
        for k, v in self.templates.iteritems():
            outpath =

    def _write_code(self):
        open(outpath, 'wt').write(template.render(subs))

    def _compile(self):
        pass

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

        from scipy.integrate import ode
        self._r = ode(self._fo_odesys.dydt, self._fo_odesys.dydt_jac)


    def integrate(self, y0, t0, tend, N, h = None, order = 0):
        y0_val_lst = [y0[k] for k in self._fo_odesys.dep_var_func_symbs]
        self._r.set_initial_value(y0_val_lst, t0)
        self._r.set_f_params(self._fo_odesys.params_val_lst)
        self._r.set_jac_params(self._fo_odesys.params_val_lst)
        if N > 0:
            # Fixed stepsize
            self._r.set_integrator('vode', method = 'bdf', with_jacobian = True)
            self.init_yout_tout_for_fixed_step_size(t0, tend, N, order)
            for i, t in enumerate(self.tout):
                self.yout[i, :] = self._r.integrate(self.tout[i])
                if order > 0:
                    self.dyout[i, :] = self.dydt(tout[i], self.yout[i, :],
                                                 self.params_val_lst)
                if order > 1:
                    self.ddyout[i,: ] = self.d2ydt2(tout[i], self.yout[i, :],
                                                 self.params_val_lst)
                assert self._r.successful()
        else:
            # Adaptive step size reporting
            # http://stackoverflow.com/questions/12926393/\
            #   using-adaptive-step-sizes-with-scipy-integrate-ode
            self._r.set_integrator('vode', method = 'bdf', with_jacobian = True, nsteps = 1)
            self._r._integrator.iwork[2] =- 1
            yout, tout= [], []
            warnings.filterwarnings("ignore", category=UserWarning)
            yout.append(y0_val_lst)
            tout.append(t0)
            while self._r.t < tend:
                self._r.integrate(tend, step=True)
                yout.append(self._r.y)
                tout.append(self._r.t)
            warnings.resetwarnings()
            self.yout = np.array(yout)
            self.tout = np.array(tout)

