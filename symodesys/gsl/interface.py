from __future__ import print_function, division

# stdlib imports
import os

# External imports
import numpy as np
import cython_gsl
from pycodeexport.codeexport import C_Code

# Intrapackage imports
from symodesys.codeexport import ODESys_Code, Binary_IVP_Integrator

class GSL_Code(ODESys_Code, C_Code):

    # Implement hash of fo_odesys and hash of code?
    # Serialization to double check against collision?

    compile_kwargs = {
        'options': ['warn', 'pic', 'fast'],
        'std': 'c99',
        'defmacros': ['GSL_RANGE_CHECK_OFF', 'HAVE_INLINE'],
        'libraries': cython_gsl.get_libraries(),
        'include_dirs': [cython_gsl.get_include(), cython_gsl.get_cython_include_dir()],
        'library_dirs': [cython_gsl.get_library_dir()]
    }

    build_files = ['prebuilt/_drivers.o',
                  'prebuilt/drivers.o',
                  'drivers.h', 'drivers.c', 'ode.h', 'Makefile', 'plot.py',
                  'symodesys_util.c', 'symodesys_util.h',
    ]

    obj_files = ['ode.o', 'drivers.o',
                 'symodesys_util.o', '_drivers.o']

    templates = ['ode_template.c',
                 'main_ex_template.c',
    ]

    source_files = ['ode.c', 'symodesys_util.c']

    so_file = '_drivers.so'

    extension_name = '_drivers'


    # ODESys_Code specific
    depv_tok = 'y' # see ode_template.c
    depv_offset = None

    param_tok = 'k' # see ode_template.c
    param_offset = None

    def __init__(self, *args, **kwargs):
        self.basedir = os.path.dirname(__file__)
        super(GSL_Code, self).__init__(*args, **kwargs)


class GSL_IVP_Integrator(Binary_IVP_Integrator):
    """
    IVP integrator using GNU Scientific Library routines odeiv2

    remember to run `clean()` on the instance on program exit. (or use IVP's context manager)
    """

    CodeClass = GSL_Code

    # step_type choices are in `step_types` in _drivers.pyx
    integrate_args = {'step_type': (
        'rk2','rk4','rkf45','rkck','rk8pd','rk2imp',
        'rk4imp','bsimp','rk1imp','msadams','msbdf'),
    }


    def _run(self, depv_init_arr, indepv_init, indepv_end, params_arr, N,
            step_type='bsimp'):
        """
        -`params`: dict with keys mathcing param_and_sol_symbs

        Notes: hmax won't be set if 0.0
        """

        assert step_type in self.integrate_args['step_type']

        if N > 0:
            # Fixed stepsize
            self.h_init = self.h_init or 1e-9
            tout, Yout = self.binary_mod.integrate_equidistant_output(
                indepv_init, indepv_end, depv_init_arr, N,
                self.h_init, self.h_max, self.abstol,
                self.reltol, params_arr, self.nderiv, step_type)
            self.tout = tout
            self.Yout = Yout
        else:
            raise NotImplementedError
