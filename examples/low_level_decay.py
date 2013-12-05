#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Example not using SimpleFirstOrderODESystem but the more low-level
FirstOrderODESystem . The main difference is that this approach
requires a bit more knowledge on how to use sympy and how symodesys
represents the odesystem internally (with great power comes great
responsibilities).
"""

# Standard library imports
from collections import OrderedDict

# External imports
import sympy
import numpy as np
import matplotlib.pyplot as plt

# Project internal imports
from symodesys import FirstOrderODESystem
from symodesys.convenience import numeric_vs_analytic
from symodesys.integrator import Mpmath_IVP_Integrator

# Global
t = sympy.symbols('t')
u = sympy.Function('u')(t)
lambda_u = sympy.symbols('lambda_u')

class Decay(FirstOrderODESystem):
    """
    Shows how to instantiate sympy.Symbol and sympy.Function instances
    manually and use them directly in FirstOrderODESystem
    """
    indepv = t
    dep_var_symbs = u,
    param_symbs   = [lambda_u]
    f = OrderedDict([(u, -lambda_u * u)])

    # Analytic sol for plotting numerical error
    @property
    def analytic_sol(self):
        from decay import analytic_decay
        return {self['u']: analytic_decay}


if __name__ == '__main__':
    numeric_vs_analytic(Decay, {'u': 1.0}, {'lambda_u': 0.2}, 0.0, 10.0, N=30,
                        integrator=Mpmath_IVP_Integrator(nderiv=2))
