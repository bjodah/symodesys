#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os

from symodesys import plot_numeric_vs_analytic
from symodesys.odepack import LSODES_IVP_Integrator
from coupled_decay import CoupledDecay

if __name__ == '__main__':
    ax = plot_numeric_vs_analytic(
        ODESys = CoupledDecay,
        y0 = {'u': 1.0, 'v': 1.0, 'w': 1.0},
        params = {'lambda_u': 1/3,'lambda_v': 1/5,'lambda_w': 1/7},
        t0 = 0,
        integrator = LSODES_IVP_Integrator(
            tempdir=os.path.join(os.path.dirname(__file__), 'build_codegen_odepack'),
            save_temp=True
            ),
        N = 100,
        tend = 5.0,
       )
