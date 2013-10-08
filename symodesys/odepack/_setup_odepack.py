#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from pycompilation import FortranCompilerRunner, pyx2obj, download_files, compile_sources

f_sources = ['opkdmain.f', 'opkda1.f', 'opkda2.f']

def main(cwd, logger):
    websrc = 'https://computation.llnl.gov/casc/odepack/software/'
    md5sums = {
        'opkda1.f':    '00a675f71ab375376bb6108d24a33c0b',
        'opkda2.f':    'dd03a71ea1a5ac746169c0279aa4c551',
        'opkdmain.f':  '47d81cc73a1e82210f47a97c43daa8cf'
    }

    download_files(websrc, f_sources, md5sums, cwd)
    # Intel Fortran fails for opkda1.f, hence prefer `gnu`
    compile_sources(f_sources, FortranCompilerRunner,
                    'prebuilt/', cwd,
                    run_linker=False,
                    options=['pic', 'warn', 'fast'],
                    preferred_vendor='gnu', metadir='prebuilt/', logger=logger)

    # Cythonize pyx file, and compile to object file
    src = 'lsodes_bdf_wrapper.pyx'
    dst = 'prebuilt/lsodes_bdf_wrapper.o'
    pyx2obj(src, dst, cwd=cwd, logger=logger, only_update=True)
