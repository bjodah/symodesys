from __future__ import print_function, division

# stdlib imports
import os

# External imports
import numpy as np
import cython_gsl
from pycodeexport.codeexport import C_Code

# Intrapackage imports
from symodesys.codeexport import ODESys_Code, Binary_IVP_Integrator

USE_LAPACK = os.environ.get('USE_LAPACK', False)
defmacros = ['SUNDIALS_DOUBLE_PRECISION']

extra_options = []
if USE_LAPACK:
    defmacros += ['USE_LAPACK']
    extra_options += ['lapack']


class CVODE_Code(ODESys_Code, C_Code):

    compile_kwargs = {
        'options': ['warn', 'pic', 'fast'] + extra_options,
        'std': 'c99',
        'defmacros': defmacros,
        'libraries': ['m', 'sundials_cvode', 'sundials_nvecserial'],
    }

    build_files = [
        'prebuilt/drivers.o',
        'prebuilt/_drivers.o',
        'drivers.h', 'drivers.c', 'Makefile', 'plot.py',
        '../shared/prebuilt/symodesys_util.o',
        'symodesys_util.h', # for main_ex
        'symodesys_util.c', # for main_ex
    ]

    obj_files = ['func.o', 'dense_jac.o', 'band_jac.o', 'drivers.o',
                 'symodesys_util.o', '_drivers.o']

    templates = ['func_template.c',
                 'dense_jac_template.c',
                 'band_jac_template.c',
                 'main_ex_template.c',
    ]

    source_files = ['func.c', 'dense_jac.c', 'band_jac.c']

    so_file = '_drivers.so'

    extension_name = '_drivers'


    # ODESys_Code specific
    depv_tok = 'NV_DATA_S(y)' # See template
    depv_offset = None

    param_tok = 'k' # See template
    param_offset = None

    def __init__(self, *args, **kwargs):
        self.basedir = os.path.dirname(__file__)
        super(CVODE_Code, self).__init__(*args, **kwargs)


class CVODE_IVP_Integrator(Binary_IVP_Integrator):
    """
    IVP integrator using GNU Scientific Library routines odeiv2

    remember to run `clean()` on the instance on program exit. (or use IVP's context manager)
    """

    CodeClass = CVODE_Code

    # step_type choices are in `step_types` in drivers_wrapper.pyx
    integrate_args = {
        'step_type': ('adams','bdf'),
        'mode': ('dense', 'band')
    }


    def _run(self, depv_init_arr, indepv_init, indepv_end, params_arr, N,
             step_type='adams', mode='dense'):
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
                self.reltol, params_arr, self.nderiv, step_type, mode)
            self.tout = tout
            self.Yout = Yout
        else:
            raise NotImplementedError
