# -*- coding: utf-8 -*-

import os
import shutil
from pycompilation import src2obj
from pycompilation.util import copy

shared_sources = ['symodesys_util.c']

def prebuild(srcdir, destdir, **kwargs):
    for cf in 'main_ex_template.c symodesys_util.h symodesys_util.c'.split():
        copy(os.path.join(srcdir, cf), destdir)
    destdir = os.path.join(destdir, 'prebuilt')
    return [src2obj(
        os.path.join(srcdir, f),
        objpath=destdir,
        options=['pic', 'warn', 'fast'],
        std='c99',
        metadir=destdir,
        only_update=True,
        **kwargs) for f in shared_sources]
