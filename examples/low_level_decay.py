#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Standard library imports
from collections import OrderedDict

# External imports
import sympy
import numpy as np
import matplotlib.pyplot as plt

# Package imports
from symodesys import FirstOrderODESystem
from symodesys.convenience import plot_numeric_error

# Global
t = sympy.symbols('t')
u = sympy.Function('u')(t)
lambda_u = sympy.symbols('lambda_u')

class Decay(FirstOrderODESystem):
    """
    Shows how to instantiate sympy.symbols and sympy.Function instances
    manually and use them directly in FirstOrderODESystem
    """
    indepv = t
    dep_var_symbs = u,
    param_symbs   = [lambda_u]
    f = OrderedDict([(u, -lambda_u * u)])

    def analytic_u(self, indep_vals, y0, params, t0):
        return y0['u'] * np.exp(-params['lambda_u']*indep_vals)

    analytic_sol = {'u': analytic_u}


if __name__ == '__main__':
    plot_numeric_error(Decay, {'u': 1.0}, {'lambda_u': 0.2}, 0.0, 10.0, N=30)
