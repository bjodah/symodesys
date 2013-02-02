#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy

import matplotlib.pyplot as plt

from symodesys.firstorder import FirstOrderODESystem
from symodesys.integrator import SciPy_IVP_Integrator

class VanDerPolOscillator(FirstOrderODESystem):

    num_dep_vars = 2
    num_params = 1

    # Following two lines are optional but useful for
    # automatic labeling when plotting:
    dep_var_symbs = sympy.symbols('u v')
    param_symbs   = sympy.symbols('mu,')

    @property
    def f(self):
        u, v = self.dep_var_symbs
        mu = self.param_symbs[0]
        return [v,
                -u + mu * v * (1 - u ** 2),
                ]


def main(mu):
    """
    """
    vdpo = VanDerPolOscillator()
    intr = SciPy_IVP_Integrator(vdpo, [mu])

    int_kwargs = {'abstol': 1e-6,
                  'reltol': 1e-6}

    #N = 0 # adaptive stepsize controls output
    N = 100
    intr.integrate([1.0, 0.0], 0.0, 10.0, N, **int_kwargs)
    intr.plot(interpolate = True)

if __name__ == '__main__':
    main(sys.argv[1])
