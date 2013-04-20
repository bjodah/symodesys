#!/usr/bin/env python
# -*- coding: utf-8 -*-

# External imports
import sympy
import numpy as np

# Package imports
from symodesys import FirstOrderODESystem
from symodesys.helpers import plot_numeric_vs_analytic

# Global
t = sympy.symbols('t')
u = sympy.Function('u')(t)
lambda_u = sympy.symbols('lambda_u')

class Decay(FirstOrderODESystem):

    dep_var_symbs = u,
    param_symbs   = lambda_u,
    f = OrderedDict([(u, -lambda_u * u)])

    def analytic_u(self, indep_vals, y0):
        return y0['u'] * np.exp(-self.param_vals_by_symb[lamba_u]*indep_vals)


if __name__ == '__main__':
    plot_numeric_vs_analytic(Decay, {'u': 1.0}, {'lambda_u': 0.2}, 10.0)
