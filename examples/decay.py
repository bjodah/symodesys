#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.firstorder import FirstOrderODESystem
from symodesys.integrator import SciPy_IVP_Integrator
from symodesys.ivp import IVP

# TODO: add use of units from Sympys physics module and enter lambda in per s, and give
#        time intervals in hours

class Decay(FirstOrderODESystem):

    #num_dep_vars = 1
    #num_params = 1

    # Following two lines are optional but useful for
    # automatic labeling when plotting:
    dep_var_tokens = 'u',
    param_tokens   = 'lambda_u',

    @property
    def f(self):
        u, = self.dep_var_func_symbs
        lambda_u, = self.param_symbs
        return {u: -lambda_u * u,
                }


def main(params):
    """
    """
    d = Decay()
    d.update_default_params_by_token(params)

    y0 = {d['u']: 1.0}
    ivp = IVP(d, y0, SciPy_IVP_Integrator)

    intgrtn_kwargs = {'abstol': 1e-6,
                      'reltol': 1e-6}

    N = 0 # adaptive stepsize controls output
    t0 = 0.0 # initial time = 0
    tend = 10.0 # final time = 10

    ivp.integrate(t0, tend, N, **intgrtn_kwargs)
    t, u = ivp.tout, ivp.yout

    plt.subplot(311)
    ivp.plot(interpolate = True, show = False)

    lambda_u = params['lambda_u']
    analytic_u = np.exp(-lambda_u*t)
    plt.subplot(312)
    plt.plot(t, (u[:, 0] - analytic_u) / intgrtn_kwargs['abstol'],
             label = 'abserr / abstol')
    plt.legend()
    plt.subplot(313)
    plt.plot(t, (u[:, 0] - analytic_u) / analytic_u / intgrtn_kwargs['reltol'],
             label = 'relerr / reltol')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        main({'lambda_u': float(sys.argv[1])})
    else:
        main({'lambda_u': 0.2})
