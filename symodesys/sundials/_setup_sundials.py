# -*- coding: utf-8 -*-

"""
Precompiles SUNDIALS wrapper
when linking don't forget libs, usually: 'm', 'sundials_cvode', 'sundials_nvecserial'
"""

import os
from pycodeexport.codeexport import make_CleverExtension_for_prebuilding_Code

def get_sundials_clever_ext(basename):
    from .interface import CVODE_Code
    return make_CleverExtension_for_prebuilding_Code(
        basename+'.sundials._drivers', CVODE_Code,
        ['drivers.c', '_drivers.pyx'],
        srcdir=os.path.join(basename, 'sundials'),
        logger=True
    )
