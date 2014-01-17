# -*- coding: utf-8 -*-

import os
from pycompilation import pyx2obj, src2obj
from pycompilation.util import copy

def prebuild(srcdir, destdir, **kwargs):
    for cf in 'Makefile drivers.c drivers.h main_ex_template.c'.split() +\
        'ode.h symodesys_util.c ode_template.c symodesys_util.h'.split():
        copy(os.path.join(srcdir, cf), destdir)
    destdir = os.path.join(destdir, 'prebuilt')
    return [
        pyx2obj(os.path.join(srcdir, 'drivers_wrapper.pyx'),
                destdir, only_update=True,
                metadir=destdir, **kwargs),
        src2obj(
            os.path.join(srcdir,'drivers.c'),
            objpath=destdir,
            defmacros=['GSL_RANGE_CHECK_OFF', 'HAVE_INLINE'],
            options=['pic', 'warn', 'fast'],
            #remember to use libs=['gsl', 'gslcblas', 'm'], when linking
            std='c99',
            metadir=destdir,
            only_update=True,
            **kwargs
        )
    ]
