# -*- coding: utf-8 -*-

import os
from pycompilation.codeexport import make_CleverExtension_for_prebuilding_Code

def get_gsl_clever_ext(basename):
    from .interface import GSL_Code
    return make_CleverExtension_for_prebuilding_Code(
        basename+'.gsl._drivers', GSL_Code,
        ['drivers.c', '_drivers.pyx'],
        srcdir=os.path.join(basename, 'gsl'),
        logger=True
    )
