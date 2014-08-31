# -*- coding: utf-8 -*-

import os
from pycodeexport.codeexport import make_PCEExtension_for_prebuilding_Code

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

def get_odepack_pce_ext(basename):
    from .interface import LSODES_Code
    return make_PCEExtension_for_prebuilding_Code(
        basename+'.odepack._drivers', LSODES_Code,
        f_sources+['_drivers.pyx'], # 'drivers.f90' uses module from ode_template
        srcdir=os.path.join(basename, 'odepack'),
        downloads=(websrc, src_md5),
        logger=True
    )
