from __future__ import print_function, division

# stdlib imports
import os

# External imports
import numpy as np
import cython_gsl

# Intrapackage imports
from symodesys.codeexport import Generic_Code, Binary_IVP_Integrator


class GSL_Code(Generic_Code):

    # Implement hash of fo_odesys and hash of code?
    # Serialization to double check against collision?

    copy_files = ['drivers.c', 'pyinterface.pyx',
                  'drivers.h', 'ode.h', 'Makefile']

    templates = ['ode_template.c',
                 'main_ex_template.c',
             ]

    _source_files = copy_files[:2] + templates[:1]

    def __init__(self, *args, **kwargs):
        self._basedir = os.path.dirname(__file__)
        super(GSL_Code, self).__init__(*args, **kwargs)
        self._include_dirs.append(cython_gsl.get_include())
        self._libraries.extend(cython_gsl.get_libraries())
        self._library_dirs.append(cython_gsl.get_library_dir())
        self._include_dirs.append(cython_gsl.get_cython_include_dir())


class GSL_IVP_Integrator(Binary_IVP_Integrator):
    """
    IVP integrator using GNU Scientific Library routines odeiv2

    remember to run `clean()` on the instance on program exit. (or use IVP's context manager)
    """

    CodeClass = GSL_Code

    # step_type choices are in `step_types` in pyinterface.pyx
    integrate_args = {'step_type': (
        'rk2','rk4','rkf45','rkck','rk8pd','rk2imp',
        'rk4imp','bsimp','rk1imp','msadams','msbdf'),
    }


    def run(self, y0, t0, tend, param_vals, N,
                  abstol=None, reltol=None, h=None, h_init=None, h_max=0.0, nderiv=None,
            **kwargs):
        """
        hmax won't be set if 0.0
        """
        for k,v in kwargs.items():
            # Assert valid option provided
            if k in self.integrate_args:
                assert v in self.integrate_args[k]
        self.nderiv = nderiv or self.nderiv
        if h_init == None:
            h_init = 1e-9 # TODO: along the lines of:
            #   h_init=calc_h_init(y0, dydt, jac, abstol, reltol)

        # Set (c-program) init vals to the possible subset of non_anlytic depv
        y0_arr = np.array(
            [y0[k] for k in self._fo_odesys.non_analytic_depv],
            dtype = np.float64)

        # Extend (c-program) params with y0 values of analytic functions
        prog_param_vals = dict(
            param_vals.items() + \
            [(k, v) for k, v in y0.items() if k \
             in self._fo_odesys.analytic_depv])

        # The C-program knows nothing about dicts, provide values as an array
        params_arr = np.array([prog_param_vals[k] for k \
                               in self._code.prog_param_symbs],
                              dtype = np.float64)
        if N > 0:
            # Fixed stepsize
            #self.init_Yout_tout_for_fixed_step_size(t0, tend, N)
            tout, Yout = self.binary.integrate_equidistant_output(
                t0, tend, y0_arr, N, h_init, h_max, self.abstol,
                self.reltol, params_arr, self.nderiv, **kwargs)
            self.tout = tout
            self.Yout = Yout
        else:
            raise NotImplementedError
