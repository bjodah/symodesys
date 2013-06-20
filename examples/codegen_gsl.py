#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os

from symodesys import plot_numeric_vs_analytic
from symodesys.gsl import GSL_IVP_Integrator
from coupled_decay import CoupledDecay

import logging



if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    ax = plot_numeric_vs_analytic(
        ODESys = CoupledDecay,
        y0 = {'u': 1.0, 'v': 1.0, 'w': 1.0},
        params = {'lambda_u': 1/3,'lambda_v': 1/5,'lambda_w': 1/7},
        t0 = 0,
        integrator = GSL_IVP_Integrator(
            tempdir=os.path.join(os.path.dirname(__file__), 'build_codegen_gsl'),
            save_temp=True,
            logger=logger,
            ),
        N = 100,
        tend = 5.0,
        show=True
       )
