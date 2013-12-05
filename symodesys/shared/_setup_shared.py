# -*- coding: utf-8 -*-

import os
import shutil
from pycompilation import src2obj

shared_sources = ['symodesys_util.c']

def main(dst, **kwargs):
    return [src2obj(
        f,
        objpath=dst,
        options=['pic', 'warn', 'fast'],
        std='c99',
        metadir=dst,
        only_update=True,
        **kwargs) for f in shared_sources]
