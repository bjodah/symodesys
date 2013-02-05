#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.firstorder import FirstOrderODESystem
from symodesys.integrator import SciPy_IVP_Integrator

# TODO
# Implement automatic resolution of N - number of chained decays via
# Bateman's equations

class CoupledDecay(FirstOrderODESystem):

    num_dep_vars = 2
    num_params = 2

    # Following two lines are optional but useful for
    # automatic labeling when plotting:
    dep_var_symbs = sympy.symbols('u v')
    param_symbs   = sympy.symbols('lambda_u, lambda_v')

    @property
    def f(self):
        u, v = self.dep_var_symbs
        lambda_u, lambda_v = self.param_symbs
        return [-lambda_u * u,
                lambda_u * u - lambda_v * v,
                ]


def main(lambda_u, lambda_v):
    """
    """
    cd = CoupledDecay()
    intr = SciPy_IVP_Integrator(cd, [lambda_u, lambda_v])

    int_kwargs = {'abstol': 1e-6,
                  'reltol': 1e-6}

    #N = 0 # adaptive stepsize controls output
    N = 100
    t0 = 0.0
    tend = 10.0
    u0 = 1.0
    v0 = 0.0

    intr.integrate([u0, v0], t0, tend, N, **int_kwargs)

    t = np.linspace(t0, tend, N)
    analytic_u = np.exp(-lambda_u*t)
    analytic_v = u0 * lambda_u / (lambda_v - lambda_u) * (np.exp(-lambda_u*t) - \
                                                       np.exp( - lambda_v * t))


    plt.subplot(311)
    intr.plot(interpolate = True, show = False)

    plt.subplot(312)
    plt.plot(intr.tout, (intr.yout[:, 0] - analytic_u) / int_kwargs['abstol'],
             label = 'u abserr / abstol')
    plt.plot(intr.tout, (intr.yout[:, 1] - analytic_v) / int_kwargs['abstol'],
             label = 'v abserr / abstol')
    plt.legend()
    plt.subplot(313)
    plt.plot(intr.tout, (intr.yout[:, 0] - analytic_u) / analytic_u / int_kwargs['reltol'],
             label = 'u relerr / reltol')
    plt.plot(intr.tout, (intr.yout[:, 1] - analytic_v) / analytic_v / int_kwargs['reltol'],
             label = 'v relerr / reltol')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main(float(sys.argv[1]), float(sys.argv[2]))
