# -*- coding: utf-8 -*-

import os
from pycompilation import pyx2obj, src2obj
from pycompilation.util import copy

USE_LAPACK = os.environ.get('USE_LAPACK', False)

def prebuild(srcdir, destdir, **kwargs):
    # Cythonize pyx file and compile to object

    defmacros = ['SUNDIALS_DOUBLE_PRECISION']
    options = ['pic', 'warn', 'fast']
    if USE_LAPACK:
        defmacros += ['USE_LAPACK']
        options += ['lapack']

    for cf in 'Makefile band_jac_template.c func_template.c'.split()+\
        'dense_jac_template.c symodesys_util.c drivers.c'.split() +\
        'symodesys_util.h drivers.h main_ex_template.c'.split():
        copy(os.path.join(srcdir, cf), destdir)

    destdir = os.path.join(destdir, 'prebuilt')
    return [
        pyx2obj(
            os.path.join(srcdir,'drivers_wrapper.pyx'),
            destdir, only_update=True, metadir=destdir, **kwargs),
        src2obj(
            os.path.join(srcdir, 'drivers.c'),
            objpath=destdir,
            options=options,
            std='c99',
            #libs=['m', 'sundials_cvode', 'sundials_nvecserial'], #when linking
            defmacros=defmacros,
            metadir=destdir,
            only_update=True,
            **kwargs)
    ]
