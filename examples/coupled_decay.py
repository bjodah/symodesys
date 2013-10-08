#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.odesys import SimpleFirstOrderODESystem
from symodesys.ivp import IVP
from symodesys.convenience import plot_numeric_error
from symodesys.integrator import SciPy_IVP_Integrator

# TODO
# Implement automatic resolution of N - number of chained decays via
# Bateman's equations <--- No, doesnt handle lambda_k = lambda_l

class CoupledDecay(SimpleFirstOrderODESystem):

    # Following two lines are optional but useful for
    # automatic labeling when plotting:
    depv_tokens = 'u v w'.split()
    param_tokens   = 'lambda_u lambda_v lambda_w'.split()

    @property
    def expressions(self):
        u, v, w = self['u'], self['v'], self['w']
        lambda_u, lambda_v, lambda_w = self.param_symbs
        return {u: -lambda_u * u,
                v: lambda_u * u - lambda_v * v,
                w: lambda_v * v - lambda_w * w,
                }

    def analytic_u(self, indep_vals, y0, params, t0):
        return y0['u'] * np.exp(-params['lambda_u']*indep_vals)


    def analytic_v(self, indep_vals, y0, params, t0):
        return y0['v'] * np.exp(-params['lambda_v'] * indep_vals) + \
                 y0['u'] * params['lambda_u'] / \
                 (params['lambda_v'] - params['lambda_u']) * \
                 (np.exp(-params['lambda_u']*indep_vals) - \
                  np.exp( - params['lambda_v'] * indep_vals))

    def analytic_w(self, indep_vals, y0, params, t0):
        return y0['w'] * np.exp(-params['lambda_w'] * indep_vals) + \
                 y0['v'] * params['lambda_v'] / \
                 (params['lambda_w'] - params['lambda_v']) * \
                 (np.exp(-params['lambda_v']*indep_vals) - \
                  np.exp(-params['lambda_w']*indep_vals)) + \
                 params['lambda_v'] * params['lambda_u'] * \
                 y0['u'] / (params['lambda_v'] - \
                            params['lambda_u']) * \
                 (1 / (params['lambda_w'] - \
                       params['lambda_u']) * \
                  (np.exp( - params['lambda_u'] * indep_vals) - \
                   np.exp( - params['lambda_w'] * indep_vals)) - \
                  1 / (params['lambda_w'] - \
                       params['lambda_v']) * \
                  (np.exp( - params['lambda_v'] * indep_vals) - \
                   np.exp( - params['lambda_w'] * indep_vals)))

    analytic_sol = {'u': analytic_u, 'v': analytic_v, 'w': analytic_w}




if __name__ == '__main__':
    if len(sys.argv) > 1:
        lambda_u = float(sys.argv[1])
        lambda_v = float(sys.argv[2])
        lambda_w = float(sys.argv[3])
    else:
        lambda_u, lambda_v, lambda_w = 3.0, 2.0, 1.0

    integrator = SciPy_IVP_Integrator()
    integrator.abstol = 1e-7
    integrator.reltol = 1e-7
    plot_numeric_error(
        CoupledDecay,
        {'u': 7.0,
         'v': 5.0,
         'w': 3.0,
     },
        params = {
            'lambda_u': lambda_u,
            'lambda_v': lambda_v,
            'lambda_w': lambda_w,
        },
        indepv_init = 0.0,
        indepv_end = 1.5,
        integrator=integrator
    )
