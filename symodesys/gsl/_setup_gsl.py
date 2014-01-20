# -*- coding: utf-8 -*-

from pycompilation.codeexport import prebuild_Code

def prebuild(srcdir, destdir, build_temp, **kwargs):
    from .interface import GSL_Code as Code

    all_sources = ['drivers.c', '_drivers.pyx']
    return prebuild_Code(
        srcdir, destdir, build_temp, Code, all_sources,
        per_file_kwargs={
            'drivers.c': {
                'defmacros': ['GSL_RANGE_CHECK_OFF', 'HAVE_INLINE'],
                'std': 'c99',
            }
        },
        **kwargs
    )
