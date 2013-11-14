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

class CVODE_Code(ODESys_Code, C_Code):

    compilation_options = ['c99'] # see pycompilation

    copy_files = ['prebuilt/drivers_wrapper.o',
                  'prebuilt/drivers.o',
                  'drivers.h', 'drivers.c', 'Makefile', 'plot.py',
                  '../shared/prebuilt/symodesys_util.o', #, 'symodesys_util.h',
                  'prebuilt/'+CCompilerRunner.metadata_filename, # <--- Ensure we use the same compiler
    ]

    obj_files = ['func.o', 'dense_jac.o', 'band_jac.o', 'drivers.o',
                 'drivers_wrapper.o', 'symodesys_util.o']

    templates = ['func_template.c',
                 'dense_jac_template.c',
                 'band_jac_template.c',
                 'main_ex_template.c',
    ]

    source_files = ['func.c', 'dense_jac.c', 'band_jac.c']

    so_file = 'drivers_wrapper.so'

    extension_name = 'drivers_wrapper'


    # ODESys_Code specific
    depv_tok = 'NV_DATA_S(y)' # See template
    depv_offset = None

    param_tok = 'k' # See template
    param_offset = None

    def __init__(self, *args, **kwargs):
        self._basedir = os.path.dirname(__file__)
        super(CVODE_Code, self).__init__(*args, **kwargs)
        # self.inc_dirs.append()
        # self.inc_dirs.append()
        self.libs.extend(['m', 'sundials_cvode', 'sundials_nvecserial'])
        self.defmacros.extend(['SUNDIALS_DOUBLE_PRECISION'])
        # self.lib_dirs.append()


class CVODE_IVP_Integrator(Binary_IVP_Integrator):
    """
    IVP integrator using GNU Scientific Library routines odeiv2

    remember to run `clean()` on the instance on program exit. (or use IVP's context manager)
    """

    CodeClass = CVODE_Code

    # step_type choices are in `step_types` in drivers_wrapper.pyx
    integrate_args = {'step_type': ('adams','bdf'),
    }


    def _run(self, depv_init_arr, indepv_init, indepv_end, params_arr, N,
            step_type='adams'):
        """
        -`params`: dict with keys mathcing param_and_sol_symbs

        Notes: hmax won't be set if 0.0
        """

        assert step_type in self.integrate_args['step_type']

        if N > 0:
            # Fixed stepsize
            self.h_init = self.h_init or 0.0 # 0.0 -> will not set -> CVODE chooses h_init
            tout, Yout = self.binary_mod.integrate_equidistant_output(
                indepv_init, indepv_end, depv_init_arr, N,
                self.h_init, self.h_max, self.abstol,
                self.reltol, params_arr, self.nderiv, step_type)
            self.tout = tout
            self.Yout = Yout
        else:
            raise NotImplementedError
