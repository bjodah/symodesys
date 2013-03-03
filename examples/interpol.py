#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.firstorder import SimpleFirstOrderODESystem, FirstOrderODESystem
from symodesys.ivp import IVP


class X5(SimpleFirstOrderODESystem):

    dep_var_tokens = 'y',
    def init_f(self):
        self.f = {self['y']: self.indep_var_symb ** 4}

    def analytic_y(self, indep_vals, t0, y0):
        return y0['y'] + indep_vals ** 5 / 5 - t0 ** 5


def plot_numeric_vs_analytic(ODESys,
                             indep_var_lim,
                             init_dep_var_vals_by_token,
                             param_vals,
                             N = 0):
    """
    Integrate
    """

    # Setup
    odesys = ODESys()
    odesys.update_params_by_token(param_vals)

    t0, tend = indep_var_lim
    y0 = {odesys[k]: v for k, v in init_dep_var_vals_by_token.items()}

    # Solve
    ivp = IVP(odesys, y0, t0)
    ivp.integrate(tend, N = N)

    # Anlyse output
    t, y = ivp.tout, ivp.yout

    plt.subplot(311)
    ivp.plot(interpolate = True, show = False)
    analytic_y = odesys.analytic_y(t, t0, init_dep_var_vals_by_token)
    plot_t = np.linspace(t0, tend, 50)
    plot_ay = odesys.analytic_y(plot_t, t0, init_dep_var_vals_by_token)
    plt.plot(plot_t, plot_ay, label = 'Analytic')

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
        ODESys = X5,
        indep_var_lim = (0, 5.0),
        init_dep_var_vals_by_token = {'y': 1.0},
        param_vals = {},
        N = 2)
