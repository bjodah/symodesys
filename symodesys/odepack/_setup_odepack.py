#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import os
import sys
import subprocess
import shutil

from pycompilation import (
    FortranCompilerRunner, CCompilerRunner, pyx2obj,
    md5_of_file, missing_or_other_newer
)

def main(cwd, logger):
    opkfiles = ['opkdmain.f', 'opkda1.f', 'opkda2.f']
    opkweb = 'https://computation.llnl.gov/casc/odepack/software/'
    md5sums = {
        'opkda1.f':    '00a675f71ab375376bb6108d24a33c0b',
        'opkda2.f':    'dd03a71ea1a5ac746169c0279aa4c551',
        'opkdmain.f':  '47d81cc73a1e82210f47a97c43daa8cf'
    }

    # Download sources ----------------------------------------
    for f in opkfiles:
        fpath = os.path.join(cwd, f)
        if not os.path.exists(fpath):
            import urllib2
            print('Downloading: {}'.format(opkweb+f))
            open(fpath, 'wt').write(urllib2.urlopen(opkweb+f).read())
        fmd5 = md5_of_file(fpath).hexdigest()
        if fmd5 != md5sums[f]:
            raise ValueError("""Warning: MD5 sum of {} differs from that provided in setup.py.
            i.e. {} vs. {}""".format(fpath, fmd5, md5sums[f]))

    # (Pre)compile sources ----------------------------------------

    # Distutils does not allow to use .o files in compilation
    # (see http://bugs.python.org/issue5372)
    # hence the compilation of ODEPACK is done once and for all and
    # saved in prebuilt dir

    for f in opkfiles:
        # compiles to: prebuilt/opkd{a1,a2,main}.o
        fpath = os.path.join(cwd, f)
        name, ext = os.path.splitext(f)
        dst = os.path.join(cwd, 'prebuilt', name+'.o') # .ext -> .o
        if missing_or_other_newer(dst, f):
            # Intel Fortran fails for opkda1.f, hence prefer `gnu`
            runner = FortranCompilerRunner(
                [fpath], dst, run_linker=False,
                cwd=cwd, options=['pic', 'warn', 'fast'],
                preferred_vendor='gnu', metadir='prebuilt/', logger=logger)
            runner.run()
        else:
            print("Found {}, did not recompile.".format(dst))


    # Cythonize pyx file
    src = 'pylsodes_bdf.pyx'
    dst = 'prebuilt/pylsodes_bdf.o'
    pyx2obj(src, dst, cwd=cwd, logger=logger, only_update=True)