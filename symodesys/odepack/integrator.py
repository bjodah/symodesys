from __future__ import print_function, division

# stdlib imports
import os

# External imports
import numpy as np
import cython_gsl

# Intrapackage imports
from symodesys.codeexport import Generic_Code, Binary_IVP_Integrator
from symodesys.helpers.compilation import FortranCompilerRunner #, CCompilerRunner

class LSODES_Code(Generic_Code):

    syntax = 'F'

    copy_files = [
        'prebuilt/opkda1.o',
        'prebuilt/opkda2.o',
        'prebuilt/opkdmain.o',
        'prebuilt/pylsodes_bdf.o',
        'types.f90',
        'lsodes_bdf.f90',
        'prebuilt/'+FortranCompilerRunner.metadata_filename, # <--- Make sure we compile with same compiler
    ]

    templates = ['lsodes_bdf_wrapper_template.f90',
                 'ode_template.f90',
        ]

    _source_files = None

    def __init__(self, *args, **kwargs):
        self._basedir = os.path.dirname(__file__)
        super(LSODES_Code, self).__init__(*args, **kwargs)
        self._cached_files = ['ode.mod', 'types.mod',
                              'lsodes_bdf.mod', 'lsodes_bdf_wrapper.mod']

    def _compile(self, extension_name='pyinterface'):
        # Compile types.o ode,o lsodes_bdf.o lsodes_bdf_wrapper.o
        for f in ['types.f90', 'ode.f90', 'lsodes_bdf.f90',
                  'lsodes_bdf_wrapper.f90']:
            runner = FortranCompilerRunner(
                f, f.replace('.f90','.o'), run_linker=False,
                cwd=self._tempdir, options=['pic', 'warn', 'fast'], verbose=True,
                preferred_vendor='gnu')
            out, err, exit_status = runner.run()
            if exit_status != 0:
                print(out)
                print(err)
            else:
                print('...Success!')

        # Generate shared object for importing:
        from distutils.sysconfig import get_config_vars
        pylibs = [x[2:] for x in get_config_vars('BLDLIBRARY')[0].split() if x.startswith('-l')]
        cc = get_config_vars('BLDSHARED')[0]
        compilername, flags = cc.split()[0], cc.split()[1:] # e.g. gcc, ['-pthread', ...
        so_file = 'pylsodes_bdf.so'
        runner = FortranCompilerRunner(
            ['opkda1.o', 'opkda2.o', 'opkdmain.o', 'types.o', 'ode.o',
             'lsodes_bdf.o', 'lsodes_bdf_wrapper.o', 'pylsodes_bdf.o'],
            so_file, flags, #,compiler=   <--- We need fortran compiler.
            cwd=self._tempdir, libs=pylibs,
            verbose=True, preferred_vendor='gnu')
        out, err, exit_status = runner.run()
        if exit_status != 0:
            print(out)
            print(err)
        else:
            print('...Success!')
        return os.path.join(self._tempdir, so_file)

class LSODES_IVP_Integrator(Binary_IVP_Integrator):
    """
    IVP integrator using the sparse LSODES solver in ODEPACK (double precision)

    remember to run `clean()` on the instance on program exit. (or use IVP's context manager)
    """

    CodeClass = LSODES_Code

    def run(self, y0, t0, tend, param_vals, N,
            abstol=None, reltol=None, h=None, h_init=None, h_max=0.0, nderiv=None):
        """
        hmax won't be set if 0.0
        """
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
                t0, tend, y0_arr, N, h_init, h_max,
                self.abstol, self.reltol, params_arr, self.nderiv)
            self.tout = tout
            self.Yout = Yout
        else:
            raise NotImplementedError
