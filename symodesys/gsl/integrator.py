from __future__ import print_function, division

# stdlib imports
import os

# External imports
import numpy as np
import cython_gsl
from pycompilation import CCompilerRunner
from pycompilation.codeexport import C_Code

# Intrapackage imports
from symodesys.codeexport import ODESys_Code, Binary_IVP_Integrator

class GSL_Code(ODESys_Code, C_Code):

    # Implement hash of fo_odesys and hash of code?
    # Serialization to double check against collision?

    _copy_files = ['prebuilt/drivers_wrapper.o',
                   'prebuilt/drivers.o',
                   'drivers.h', 'drivers.c', 'ode.h', 'Makefile', 'plot.py',
                   'prebuilt/'+CCompilerRunner.metadata_filename, # <--- Make sure we compile with same compiler
               ]

    _obj_files = ['ode.o', 'drivers.o', 'drivers_wrapper.o']

    _templates = ['ode_template.c',
                 'main_ex_template.c',
              ]

    _source_files = ['ode.c']

    _so_file = 'drivers_wrapper.so'

    extension_name = 'drivers_wrapper'


    # ODESys_Code specific
    depv_tok = 'y' # see ode_template.c
    depv_offset = None

    param_tok = 'k' # see ode_template.c
    param_offset = None

    def __init__(self, *args, **kwargs):
        self._basedir = os.path.dirname(__file__)
        super(GSL_Code, self).__init__(*args, **kwargs)
        self._include_dirs.append(cython_gsl.get_include())
        self._include_dirs.append(cython_gsl.get_cython_include_dir())
        self._libraries.extend(cython_gsl.get_libraries())
        self._library_dirs.append(cython_gsl.get_library_dir())


class GSL_IVP_Integrator(Binary_IVP_Integrator):
    """
    IVP integrator using GNU Scientific Library routines odeiv2

    remember to run `clean()` on the instance on program exit. (or use IVP's context manager)
    """

    CodeClass = GSL_Code

    # step_type choices are in `step_types` in drivers_wrapper.pyx
    integrate_args = {'step_type': (
        'rk2','rk4','rkf45','rkck','rk8pd','rk2imp',
        'rk4imp','bsimp','rk1imp','msadams','msbdf'),
    }


    def run(self, y0, t0, tend, params, N,
            abstol=None, reltol=None, h=None,
            h_init=None, h_max=0.0, nderiv=None,
            step_type='bsimp'):
        """
        -`params`: dict with keys mathcing param_and_sol_symbs

        Notes: hmax won't be set if 0.0
        """
        # The C-program knows nothing about dicts, provide values as an array
        # Set (c-program) init vals to the (subset of non_anlytic) depv
        y0_arr = np.array(
            [y0[k] for k in self._fo_odesys.na_depv],
            dtype = np.float64)

        params_arr = np.array([params[k] for k \
                               in self._code.prog_param_symbs],
                              dtype = np.float64)

        assert step_type in self.integrate_args['step_type']

        # Below is a somewhat ad-hoc sanity check of condition
        # of jacobian in the starting point, even though some
        # solvers might not need the jacobian at the starting point.
        # To work around this one would ideally use a variable transformation
        # and/or solving/estimating parts of the problem analytically
        jac_cond = np.linalg.cond(self._fo_odesys.evaluate_na_jac(
            t0, y0_arr, params_arr))
        print(jac_cond)
        if jac_cond*np.finfo(np.float64).eps > max(self.abstol, self.reltol):
            raise RuntimeError(("Unlikely that Jacboian with condition: {} "+\
                               "will work with requested tolerances.").format(
                                   jac_cond))

        self.nderiv = nderiv or self.nderiv
        if h_init == None:
            h_init = 1e-9 # TODO: along the lines of:
            #   h_init=calc_h_init(y0, dydt, jac, abstol, reltol)

        if N > 0:
            # Fixed stepsize
            #self.init_Yout_tout_for_fixed_step_size(t0, tend, N)
            tout, Yout = self.binary_mod.integrate_equidistant_output(
                t0, tend, y0_arr, N, h_init, h_max, self.abstol,
                self.reltol, params_arr, self.nderiv, step_type)
            self.tout = tout
            self.Yout = Yout
        else:
            raise NotImplementedError
