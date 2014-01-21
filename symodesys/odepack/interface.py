from __future__ import print_function, division

# stdlib imports
import os

# External imports
import numpy as np
import cython_gsl
from pycompilation import FortranCompilerRunner
from pycompilation.codeexport import F90_Code
# Intrapackage imports
from symodesys.codeexport import ODESys_Code, Binary_IVP_Integrator

class LSODES_Code(ODESys_Code, F90_Code):

    compilation_options = {
        'options': ['warn', 'pic', 'fast'],
    }

    # the files from prebuilt/ are built by _setup_odepack.py
    build_files = [
        'prebuilt/opkda1.o',
        'prebuilt/opkda2.o',
        'prebuilt/opkdmain.o',
        'prebuilt/_drivers.o',
        'drivers.f90',
        'prebuilt/'+FortranCompilerRunner.metadata_filename, # <-- Ensure same compiler
    ]

    templates = ['ode_template.f90']

    source_files = ['ode.f90', 'drivers.f90']

    obj_files = ['opkda1.o', 'opkda2.o', 'opkdmain.o',
                  'ode.o', 'drivers.o',
                  '_drivers.o']

    so_file = '_drivers.so'

    preferred_vendor = 'gnu'


    # ODESys_Code specific

    depv_tok = 'y' # see ode_template.f90
    depv_offset = 1 # fortran arrays starts at 1
    # Oddity of ODEPACK, params passed at end of y-array:
    param_tok = 'y' # see ode_template.f90
    param_offset = property(lambda self: 1+self.NY)

    def __init__(self, *args, **kwargs):
        self._basedir = self._basedir or os.path.dirname(__file__)
        super(LSODES_Code, self).__init__(*args, **kwargs)


class LSODES_IVP_Integrator(Binary_IVP_Integrator):
    """
    IVP integrator using the sparse LSODES solver in ODEPACK (double precision)

    remember to run `clean()` on the instance on program exit. (or use IVP's context manager)
    """

    CodeClass = LSODES_Code

    def _run(self, depv_init_arr, indepv_init, indepv_end, params_arr, N):
        if N > 0:
            # Fixed stepsize
            self.tout, self.Yout = self.binary_mod.integrate_equidistant_output(
                indepv_init, indepv_end, depv_init_arr, N, self.h_init, self.h_max,
                self.abstol, self.reltol, params_arr, self.nderiv)
        else:
            raise NotImplementedError
