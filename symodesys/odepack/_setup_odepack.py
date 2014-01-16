#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import os

from pycompilation import pyx2obj, download_files, compile_sources
from pycompilation.util import copy

from .interface import LSODES_Code


websrc = 'https://computation.llnl.gov/casc/odepack/software/'
src_md5 = {
    'opkda1.f':    '00a675f71ab375376bb6108d24a33c0b',
    'opkda2.f':    'dd03a71ea1a5ac746169c0279aa4c551',
    'opkdmain.f':  '47d81cc73a1e82210f47a97c43daa8cf'
}

f_sources = src_md5.keys()

def prebuild(srcdir, destdir, **kwargs):

    download_files(websrc, f_sources, src_md5,
                   cwd=srcdir)
    for cf in filter(lambda x: not x.startswith('prebuilt'),
                     LSODES_Code.copy_files):
        copy(os.path.join(srcdir, cf), destdir)
    destdir = os.path.join(destdir, 'prebuilt')
    return compile_sources(
        map(lambda x: os.path.join(srcdir, x), f_sources),
        destdir=destdir,
        run_linker=False,
        options=['pic', 'warn', 'fast'],
        preferred_vendor='gnu', # ifort chokes on opkda1.f
        only_update=True,
        metadir=destdir,
        **kwargs) + [pyx2obj(
            os.path.join(srcdir, 'lsodes_bdf_wrapper.pyx'),
            destdir, only_update=True, metadir=destdir, **kwargs)]
