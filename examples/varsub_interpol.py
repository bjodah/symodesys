#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from decay import Decay

from interpol import X5

# TODO: Let IVP do transfm and back-transfm internally + convert PiecewisePolynomial


def plot_numeric_vs_analytic(Sys, indep_var_lim,
                             init_dep_var_vals_by_token, param_vals, N = 0):
    """
    Integrate
    """

    # Setup
    odesys = ODESys()
    param_vals_by_symb = odesys.get_param_vals_by_symb_from_by_token(
        param_vals)

    t0, tend = indep_var_lim
    y0 = {odesys[k]: v for k, v in init_dep_var_vals_by_token.items()}

    # Instantiate Initial Value Problem:
    ivp = IVP(odesys, y0, param_vals_by_symb, t0)

    # Set IVP to interally use variable substitution:
    # z=log(u)
    z = sympy.Function('z')(odesys.indepv)
    u = odesys['u']
    trnsfm = {z: sympy.log(u)}
    inv_trsfm = {u: sympy.exp(z)}
    # This call also converts initial values:
    ivp.use_internal_trnsfm(trnsfm, inv_trsfm)
    y0 = {odesys[dv]: y0[dv] if dv in y0 else trnsfm[dv].subs(y0) for dv in odesys.all_depv}

    # Solve IVP
    ivp.integrate(tend, N = N)

    # Anlyse output
    t, y = ivp.tout, ivp.yout

    ivp.plot(interpolate = True, show = False)
    analytic_y = odesys.analytic_y(t, t0, init_dep_var_vals_by_token)
    plot_t = np.linspace(t0, tend, 50)
    plot_ay = odesys.analytic_y(plot_t, t0, init_dep_var_vals_by_token)
    plt.plot(plot_t, plot_ay, ls = '--', label = 'Analytic')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    plot_numeric_vs_analytic(
        Sys = Decay,
        indep_var_lim = (0, 10.0),
        init_dep_var_vals_by_token = {'u': 1.0},
        param_vals = {'lambda_u': 0.2},
        N = 0)
