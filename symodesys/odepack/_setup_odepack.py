#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

from pycompilation import pyx2obj, download_files, compile_sources


websrc = 'https://computation.llnl.gov/casc/odepack/software/'
src_md5 = {
    'opkda1.f':    '00a675f71ab375376bb6108d24a33c0b',
    'opkda2.f':    'dd03a71ea1a5ac746169c0279aa4c551',
    'opkdmain.f':  '47d81cc73a1e82210f47a97c43daa8cf'
}

f_sources = src_md5.keys()

def main(dst, **kwargs):

    download_files(websrc, f_sources, src_md5,
                   cwd=kwargs.get('cwd','.'))
    return compile_sources(
        f_sources, destdir=dst,
                    run_linker=False,
        options=['pic', 'warn', 'fast'],
        preferred_vendor='gnu', # ifort chokes on opkda1.f
        only_update=True,
        metadir=dst,
        **kwargs) + [pyx2obj('lsodes_bdf_wrapper.pyx', dst,
                             only_update=True, metadir=dst, **kwargs)]
