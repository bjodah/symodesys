#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Demonstrates codeexport, compilation and
integration using LSODES from ODEPACK
"""

from __future__ import division

import os
import logging

from symodesys.odepack import LSODES_IVP_Integrator

from coupled_decay import coupled_decay_numeric_vs_analytic

if __name__ == '__main__':
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    coupled_decay_numeric_vs_analytic(
        LSODES_IVP_Integrator,
        {
            'tempdir': os.path.join(
                os.path.dirname(__file__),
                'build_codegen_odepack'
            ),
            'save_temp': True,
            'logger': logger,
        },
        logger=logger,
    )
