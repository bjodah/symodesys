# -*- coding: utf-8 -*-

from pycompilation.codeexport import prebuild_Code

"""
Precompiles ODEPACK sources from netlib (downloaded when needed)
to object files for speeding up compilation/linking on the fly
during runtime.
"""

websrc = 'https://computation.llnl.gov/casc/odepack/software/'
src_md5 = {
    'opkda1.f':    '00a675f71ab375376bb6108d24a33c0b',
    'opkda2.f':    'dd03a71ea1a5ac746169c0279aa4c551',
    'opkdmain.f':  '47d81cc73a1e82210f47a97c43daa8cf'
}


f_sources = src_md5.keys()

def prebuild(srcdir, destdir, build_temp, **kwargs):
    from .interface import LSODES_Code as Code
    all_sources = f_sources+['_lsodes_bdf.pyx']
    return prebuild_Code(
        srcdir, destdir, build_temp, Code, all_sources,
        downloads=(websrc, src_md5),
        preferred_vendor='gnu', # ifort chokes on opkda1.f
        **kwargs
    )
