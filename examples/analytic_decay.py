#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from decay import Decay




def plot_numeric_vs_analytic(Sys, indep_var_lim,
                             init_dep_var_vals_by_token, param_vals, N = 0):
    """
    Integrate
    """
    sys = Sys()
    sys.update_params_by_token(param_vals)

    y0 = {sys[k]: v for k, v in init_dep_var_vals_by_token.items()}
    print y0
    ivp = IVP(sys, y0)

    # Attempt analytic reduction
    print ivp.recursive_analytic_reduction()

    t0, tend = indep_var_lim

    ivp.integrate(t0, tend, N)
    t, y = ivp.tout, ivp.yout

    plt.subplot(311)
    ivp.plot(interpolate = True, show = False)

    analytic_y = sys.analytic_y(t, init_dep_var_vals_by_token)
    plt.subplot(312)
    plt.plot(t, y[:, 0] - analytic_y,
             label = 'abserr')
    plt.legend()
    plt.subplot(313)
    plt.plot(t, (y[:, 0] - analytic_y) / analytic_y,
             label = 'relerr')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    plot_numeric_vs_analytic(
        Sys = Decay,
        indep_var_lim = (0, 10.0),
        init_dep_var_vals_by_token = {'u': 1.0},
        param_vals = {'lambda_u': 0.2},
        N = 0)
