#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Demonstrates codeexport, compilation and
integration using Gnu Scientific Library (GSL)
"""

from __future__ import division

import os
import logging

from symodesys.gsl import GSL_IVP_Integrator

from coupled_decay import coupled_decay_numeric_vs_analytic

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    integrator_kwargs = {
        'tempdir': os.path.join(
            os.path.dirname(__file__),
            'build_codegen_gsl'
        ),
        'save_temp': True,
        'logger': logger,
    }
    coupled_decay_numeric_vs_analytic(
        GSL_IVP_Integrator,
        integrator_kwargs,
        logger=logger
    )
    fmtstr = 'To look at the code (and compile Python independent executable) look in {}'
    print(fmtstr.format(integrator_kwargs['tempdir']))
