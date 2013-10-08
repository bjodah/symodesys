#!/usr/bin/env python
# -*- coding: utf-8 -*-

# External imports
import numpy as np
from sympy import sin
from sympy import exp as e


# Package imports
from symodesys import SimpleFirstOrderODESystem
from symodesys.convenience import plot_numeric_error

"""
Product of early morning flight, not a single one is working
properly. Need to fix this when not half asleep
"""

class Sin(SimpleFirstOrderODESystem):
    depv_tokens = 'u',

    @property
    def expressions(self):
        u = self['u']
        return {u: sin(u)}

    def analytic_u(self, indep_vals, y0, params, t0):
        return -np.cos(indep_vals) + y0['u']

    analytic_sol = {'u': analytic_u}


class Nonlinear(SimpleFirstOrderODESystem):
    """
    From Kiusaalas p. 255 ex 7.5, not working correctly atm.
    """
    depv_tokens = 'u',
    param_tokens   = 'lambda_u',

    @property
    def expressions(self):
        u, l = self['u'], self['lambda_u']
        return {u: 3*u-4*e(-u)}

    def analytic_u(self, indep_vals, y0, params, t0):
        return (y0['u']-1) * np.exp(params['lambda_u']*indep_vals) + np.exp(-indep_vals)

    analytic_sol = {'u': analytic_u}

class Nonlinear2(SimpleFirstOrderODESystem):
    """
    From Kiusaalas p. 248 ex 7.2, not working correctly atm.
    """
    depv_tokens = 'u', 'up'


    @property
    def expressions(self):
        u, up, t = self['u'], self['up'], self.indepv
        return {u: up, up: -0.1*u-t}

    def analytic_u(self, indep_vals, y0, params, t0):
        return 100*indep_vals - 5*indep_vals**2 + 990*(np.exp(-0.1*indep_vals)-1)

    analytic_sol = {'u': analytic_u}


if __name__ == '__main__':
    plot_numeric_error(Nonlinear2, {'u': 0.0, 'up': 1.0}, {}, 0.0, 10.0, N=30)
    plot_numeric_error(Sin, {'u': 1.0}, {}, 0.0, 10.0, N=30)
    plot_numeric_error(Nonlinear, {'u': 1.0}, {'lambda_u': 0.2}, 0.0, 10.0, N=30)
