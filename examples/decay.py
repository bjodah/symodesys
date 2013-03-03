#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.firstorder import SimpleFirstOrderODESystem, FirstOrderODESystem
#from symodesys.integrator import SciPy_IVP_Integrator
from symodesys.ivp import IVP

# TODO: add use of units from Sympys physics module and enter lambda in per s, and give
#        time intervals in hours

u = sympy.symbols('u')
lambda_u = sympy.symbols('lambda_u')

class LowLevelDecay(FirstOrderODESystem):

    dep_var_symbs = u,
    param_symbs   = lambda_u,
    f = {u: -lambda_u * u}

    def analytic_y(self, indep_vals, y0):
        return y0['u'] * np.exp(-self.param_vals_by_symb[lamba_u]*indep_vals)


class Decay(SimpleFirstOrderODESystem):

    dep_var_tokens = 'u',
    param_tokens   = 'lambda_u',

    def init_f(self):
        self.f = {self['u']: -self['lambda_u'] * self['u']}

    def analytic_y(self, indep_vals, y0):
        return y0['u'] * np.exp(-self.params_by_token['lambda_u']*indep_vals)


def plot_numeric_vs_analytic(Sys, indep_var_lim,
                             init_dep_var_vals_by_token, param_vals, N = 0):
    """
    Integrate
    """
    sys = Sys()
    sys.update_params_by_token(param_vals)

    y0 = {sys[k]: v for k, v in init_dep_var_vals_by_token.items()}
    t0, tend = indep_var_lim
    ivp = IVP(sys, y0, t0)

    ivp.integrate(tend, N = N)
    t, y = ivp.tout, ivp.yout

    plt.subplot(311)
    ivp.plot(interpolate = True, show = False)

    analytic_y = sys.analytic_y(t, init_dep_var_vals_by_token)
    plt.subplot(312)
    plt.plot(t, (y[:, 0] - analytic_y) / ivp._integrator.abstol,
             label = 'abserr / abstol')
    plt.legend()
    plt.subplot(313)
    plt.plot(t, (y[:, 0] - analytic_y) / analytic_y / ivp._integrator.reltol,
             label = 'relerr / reltol')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    plot_numeric_vs_analytic(
        Sys = Decay,
        indep_var_lim = (0, 10.0),
        init_dep_var_vals_by_token = {'u': 1.0},
        param_vals = {'lambda_u': 0.2},
        N = 0)
