#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.firstorder import FirstOrderODESystem
from symodesys.integrator import SciPy_IVP_Integrator

class Decay(FirstOrderODESystem):

    num_dep_vars = 1
    num_params = 1

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
    d = Decay(params)
    intr = SciPy_IVP_Integrator(d)

    int_kwargs = {'abstol': 1e-6,
                  'reltol': 1e-6}

    N = 0 # adaptive stepsize controls output
    t0 = 0.0
    tend = 10.0
    y0 = {d['u']: 1.0}

    intr.integrate(y0, t0, tend, N, **int_kwargs)

    plt.subplot(311)
    intr.plot(interpolate = True, show = False)

    lambda_u = params['lambda_u']
    analytic_u = np.exp(-lambda_u*intr.tout)
    plt.subplot(312)
    plt.plot(intr.tout, (intr.yout[:, 0] - analytic_u) / int_kwargs['abstol'],
             label = 'abserr / abstol')
    plt.legend()
    plt.subplot(313)
    plt.plot(intr.tout, (intr.yout[:, 0] - analytic_u) / analytic_u / int_kwargs['reltol'],
             label = 'relerr / reltol')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        main({'lambda_u': float(sys.argv[1])})
    else:
        main({'lambda_u': 0.2})
