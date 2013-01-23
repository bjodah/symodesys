#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

from symodesys.firstorder import FirstOrderODESystem
from symodesys.integrator import IVP_ODE_Integrator

class VanDerPolOscillator(FirstOrderODESystem):

    num_dep_vars = 2
    num_params = 1

    @property
    def f(self):
        u, v = *self.dep_var_symbs
        mu = *self.param_symbs
        return [v,
                -u + mu * v * (1 - u ** 2),
                ]


vdpo = VanDerPolOscillator()

intr = IVP_ODE_Integrator(vdpo)

tols = {'abstol': 1e-8,
        'reltol': 1e-8}

tout, yout = intr.integrate([1.0, 0.0], 0.0, 100.0, 0, **tols)

plt.plot(t, y)
plt.show()
