#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os

from symodesys.convenience import assert_acceptable_numeric_errors
from symodesys.gsl import GSL_IVP_Integrator

from coupled_decay import CoupledDecay

import logging

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__name__)

    assert_acceptable_numeric_errors(
        ODESys=CoupledDecay,
        depv_init={'u': 1.0, 'v': 1.0, 'w': 1.0},
        params={'lambda_u': 1/3,'lambda_v': 1/5,'lambda_w': 1/7},
        indepv_init=0,
        integrator=GSL_IVP_Integrator(
            tempdir=os.path.join(os.path.dirname(__file__), 'build_codegen_gsl'),
            save_temp=True,
            logger=logger,
            ),
        N=100,
        indepv_end=1.0,
       )
