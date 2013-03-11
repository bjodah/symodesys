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
    odesys = Sys()
    param_vals_by_symb = odesys.get_param_vals_by_symb_from_by_token(
        param_vals)

    y0 = {odesys[k]: v for k, v in init_dep_var_vals_by_token.items()}
    t0, tend = indep_var_lim
    ivp = IVP(odesys, y0, param_vals_by_symb, t0)

    # Attempt analytic reduction
    print ivp.recursive_analytic_reduction()

    ivp.integrate(tend, N = N)
    t, u = ivp.tout, ivp.yout

    plt.subplot(311)
    ivp.plot()

    analytic_u = odesys.analytic_u(
        t, init_dep_var_vals_by_token, param_vals_by_symb)
    plt.subplot(312)
    plt.plot(t, u[:, 0] - analytic_u,
             label = 'abserr')
    plt.legend()
    plt.subplot(313)
    plt.plot(t, (u[:, 0] - analytic_u) / analytic_u,
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
