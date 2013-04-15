#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os

from symodesys import plot_numeric_vs_analytic
from symodesys.gsl import GSL_IVP_Integrator
from coupled_decay import CoupledDecay

if __name__ == '__main__':
    plot_numeric_vs_analytic(
        ODESys = CoupledDecay,
        t0 = 0,
        tend = 5.0,
        y0 = {'u': 1.0, 'v': 1.0, 'w': 1.0},
        params = {'lambda_u': 1/3,'lambda_v': 1/5,'lambda_w': 1/7},
        N = 100,
        ivp_kwargs = {'Integrator':GSL_IVP_Integrator,
                      'tempdir': os.path.join(os.path.dirname(__file__), 'codegen_out'),
                      'save_temp': True
                      }
        )
