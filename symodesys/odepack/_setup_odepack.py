#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from pycompilation import FortranCompilerRunner,  pyx2obj
from pycompilation.helpers import download_files

def main(cwd, logger):
    websrc = 'https://computation.llnl.gov/casc/odepack/software/'
    files = ['opkdmain.f', 'opkda1.f', 'opkda2.f']
    md5sums = {
        'opkda1.f':    '00a675f71ab375376bb6108d24a33c0b',
        'opkda2.f':    'dd03a71ea1a5ac746169c0279aa4c551',
        'opkdmain.f':  '47d81cc73a1e82210f47a97c43daa8cf'
    }

    download_files(websrc, files, cwd, md5sums)
    # Intel Fortran fails for opkda1.f, hence prefer `gnu`
    compile_sources(FortranCompilerRunner, files, cwd,
                    'prebuilt/', run_linker=False,
                    cwd=cwd, options=['pic', 'warn', 'fast'],
                    preferred_vendor='gnu', metadir=dst, logger=logger)

    # Cythonize pyx file, and compile to object file
    src = 'pylsodes_bdf.pyx'
    dst = 'prebuilt/pylsodes_bdf.o'
    pyx2obj(src, dst, cwd=cwd, logger=logger, only_update=True)
