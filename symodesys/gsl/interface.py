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


    def _run(self, depv_init_arr, indepv_init, indepv_end, params_arr, N,
            step_type='bsimp'):
        """
        -`params`: dict with keys mathcing param_and_sol_symbs

        Notes: hmax won't be set if 0.0
        """

        assert step_type in self.integrate_args['step_type']

        if N > 0:
            # Fixed stepsize
            #self.init_Yout_tout_for_fixed_step_size(t0, tend, N)
            tout, Yout = self.binary_mod.integrate_equidistant_output(
                indepv_init, indepv_end, depv_init_arr, N,
                self.h_init, self.h_max, self.abstol,
                self.reltol, params_arr, self.nderiv, step_type)
            self.tout = tout
            self.Yout = Yout
        else:
            raise NotImplementedError
