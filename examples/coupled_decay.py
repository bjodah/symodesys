#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy

import matplotlib.pyplot as plt

from symodesys.firstorder import FirstOrderODESystem
from symodesys.integrator import SciPy_IVP_Integrator

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


def main(lambda_u, labmda_v):
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
    v0 = 1.0

    intr.integrate([u0, v0], t0, tend, N, **int_kwargs)
    intr.plot(interpolate = True)

    t = np.linspace(t0, tend, N)
    analytic_u = np.exp(-k*t)
    #analytic_v = np.exp(-k*t)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
