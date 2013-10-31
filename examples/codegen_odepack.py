#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os, logging

from symodesys.convenience import plot_numeric_error
from symodesys.odepack import LSODES_IVP_Integrator
from coupled_decay import CoupledDecay

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    ax = plot_numeric_error(
        ODESys = CoupledDecay,
        depv_init = {'u': 1.0, 'v': 1.0, 'w': 1.0},
        params = {'lambda_u': 1/3,'lambda_v': 1/5,'lambda_w': 1/7},
        indepv_init = 0,
        integrator = LSODES_IVP_Integrator(
            tempdir=os.path.join(os.path.dirname(__file__),
                                 'build_codegen_odepack'),
            save_temp=True,
            logger=logger,
            ),
        N = 100,
        indepv_end = 5.0
       )
