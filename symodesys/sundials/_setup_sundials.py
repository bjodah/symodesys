# -*- coding: utf-8 -*-

from pycompilation.codeexport import prebuild_Code
import os

"""
Precompiles SUNDIALS wrapper
when linking don't forget libs, usually: 'm', 'sundials_cvode', 'sundials_nvecserial'
"""

USE_LAPACK = os.environ.get('USE_LAPACK', False)

def prebuild(srcdir, destdir, build_temp, **kwargs):
    defmacros = ['SUNDIALS_DOUBLE_PRECISION', 'USE_LAPACK']
    if USE_LAPACK:
        defmacros += ['USE_LAPACK']
        options += ['lapack']

    from .interface import CVODE_Code as Code
    all_sources = ['drivers.c', '_drivers.pyx']
    return prebuild_Code(
        srcdir, destdir, build_temp, Code, all_sources,
        per_file_kwargs={
            'drivers.c': {
                'defmacros': defmacros,
                'std': 'c99',
            }
        },
        **kwargs
    )
