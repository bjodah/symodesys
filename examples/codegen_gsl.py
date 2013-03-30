#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.odesys import SimpleFirstOrderODESystem, FirstOrderODESystem
from symodesys.ivp import IVP
from symodesys.gsl import GSL_IVP_Integrator


from interpol import X5



def plot_numeric_vs_analytic(ODESys,
                             indep_var_lim,
                             init_dep_var_vals_by_token,
                             param_vals,
                             N = 0):
    """
    Run compiled integrator (pointless for such a small ODE but
    useful for very large systems)
    """

    # Setup
    odesys = ODESys()
    param_vals_by_symb = odesys.get_param_vals_by_symb_from_by_token(
        param_vals)

    y0 = {odesys[k]: v for k, v in init_dep_var_vals_by_token.items()}
    t0, tend = indep_var_lim

    # Solve
    ivp = IVP(odesys, y0, param_vals_by_symb, t0,
              Integrator = GSL_IVP_Integrator,
              tempdir = os.path.join(os.path.dirname(__file__),
                                     'codegen_out'),
              save_temp = True)
    ivp.integrate(tend, N = N)

    # Anlyse output
    t, y = ivp.tout, ivp.yout

    analytic_y = odesys.analytic_y(t, t0, init_dep_var_vals_by_token)
    plot_t = np.linspace(t0, tend, 50)
    plot_ay = odesys.analytic_y(plot_t, t0, init_dep_var_vals_by_token)

    ivp.plot(interpolate = True, show = False)
    plt.plot(plot_t, plot_ay, label = 'Analytic')
    plt.show()


if __name__ == '__main__':
    plot_numeric_vs_analytic(
        ODESys = X5,
        indep_var_lim = (0, 5.0),
        init_dep_var_vals_by_token = {'y': 1.0},
        param_vals = {},
        N = 2)
