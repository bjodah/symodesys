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
    dep_var_symbs = sympy.symbols('u'),
    param_symbs   = sympy.symbols('lambda_u'),

    @property
    def f(self):
        u, = self.dep_var_symbs
        lambda_u, = self.param_symbs
        return [-lambda_u * u,
                ]


def main(lambda_u):
    """
    """
    d = Decay()
    intr = SciPy_IVP_Integrator(d, [lambda_u])

    int_kwargs = {'abstol': 1e-6,
                  'reltol': 1e-6}

    N = 0 # adaptive stepsize controls output
    t0 = 0.0
    tend = 10.0
    u0 = 1.0

    intr.integrate([u0], t0, tend, N, **int_kwargs)

    plt.subplot(311)
    intr.plot(interpolate = True, show = False)

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
    main(float(sys.argv[1]))
