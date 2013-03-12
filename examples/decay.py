#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.odesys import FirstOrderODESystem, SimpleFirstOrderODESystem
from symodesys.ivp import IVP


# TODO: add use of units from Sympys physics module and enter lambda in per s, and give
#        time intervals in hours

# TODO: Determine wheter to use:
#u = sympy.symbols('u')
#    or
#u = sympy.Function('u')(indepv)

u = sympy.symbols('u')
lambda_u = sympy.symbols('lambda_u')

class LowLevelDecay(FirstOrderODESystem):

    dep_var_symbs = u,
    param_symbs   = lambda_u,
    f = {u: -lambda_u * u}

    def analytic_u(self, indep_vals, y0):
        return y0['u'] * np.exp(-self.param_vals_by_symb[lamba_u]*indep_vals)


class Decay(SimpleFirstOrderODESystem):

    dep_var_tokens = 'u',
    param_tokens   = 'lambda_u',

    def init_f(self):
        self.f = {self['u']: -self['lambda_u'] * self['u']}

    def analytic_u(self, indep_vals, y0, param_vals):
        return y0['u'] * np.exp(-param_vals[self['lambda_u']]*indep_vals)


def plot_numeric_vs_analytic(Sys, indep_var_lim,
                             init_dep_var_vals_by_token,
                             param_vals, N = 0):
    """
    Integrate
    """
    odesys = Sys()
    param_vals_by_symb = odesys.get_param_vals_by_symb_from_by_token(
        param_vals)

    y0 = {odesys[k]: v for k, v in init_dep_var_vals_by_token.items()}
    t0, tend = indep_var_lim
    ivp = IVP(odesys, y0, param_vals_by_symb, t0)

    ivp.integrate(tend, N = N)
    t, y = ivp.tout, ivp.yout

    plt.subplot(311)
    ivp.plot(interpolate = True, show = False)

    analytic_u = odesys.analytic_u(t, init_dep_var_vals_by_token,
                                param_vals_by_symb)
    plt.subplot(312)
    plt.plot(t, (y[:, 0] - analytic_u) / ivp._integrator.abstol,
             label = 'abserr / abstol')
    plt.legend()
    plt.subplot(313)
    plt.plot(t, (y[:, 0] - analytic_u) / analytic_u / ivp._integrator.reltol,
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
