#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.firstorder import FirstOrderODESystem
from symodesys.ivp import IVP

from coupled_decay import CoupledDecay

# TODO
# Verify efficiency of general method compared to special use of
# Bateman's equations


def main(params_by_token):
    """
    """
    cd = CoupledDecay()
    param_vals_by_symb = cd.get_param_vals_by_symb_from_by_token(
        params_by_token)

    u, v, w = cd.dep_var_func_symbs

    N = 100
    t0 = 0.0
    tend = 1.5
    u0 = 7.0
    v0 = 5.0
    w0 = 3.0

    y0 = {cd['u']: u0,
          cd['v']: v0,
          cd['w']: w0,
          }

    print y0
    ivp = IVP(cd, y0, param_vals_by_symb, t0)
    new_params = ivp.recursive_analytic_reduction()
    print new_params

    ivp._Integrator.abstol = 1e-6
    ivp._Integrator.reltol = 1e-6

    ivp.integrate(tend, N = N)

    t = ivp.tout
    print y0
    analytic_u = cd.analytic_u(t, y0, param_vals_by_symb)
    analytic_v = cd.analytic_v(t, y0, param_vals_by_symb)
    analytic_w = cd.analytic_w(t, y0, param_vals_by_symb)

    uout = ivp.get_yout_by_symb(u)
    vout = ivp.get_yout_by_symb(v)
    wout = ivp.get_yout_by_symb(w)

    plt.subplot(311)
    plt.plot(t, uout, '*', label = 'Numerical u')
    plt.plot(t, vout, 'o', label = 'Numerical v')
    plt.plot(t, wout, 'd', label = 'Numerical w')
    plt.plot(t, analytic_u, label = 'Analytic u')
    plt.plot(t, analytic_v, label = 'Analytic v')
    plt.plot(t, analytic_w, label = 'Analytic w')
    plt.legend()

    plt.subplot(312)
    plt.plot(t, (uout - analytic_u) / ivp._Integrator.abstol,
             label = 'u abserr / abstol')
    plt.plot(t, (vout - analytic_v) / ivp._Integrator.abstol,
             label = 'v abserr / abstol')
    plt.plot(t, (wout - analytic_w) / ivp._Integrator.abstol,
             label = 'w abserr / abstol')
    plt.legend()

    plt.subplot(313)
    plt.plot(t, (uout - analytic_u) / analytic_u / ivp._Integrator.reltol,
             label = 'u relerr / reltol')
    plt.plot(t, (vout - analytic_v) / analytic_v / ivp._Integrator.reltol,
             label = 'v relerr / reltol')
    plt.plot(t, (wout - analytic_w) / analytic_w / ivp._Integrator.reltol,
             label = 'w relerr / reltol')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        lambda_u, lambda_v, lambda_w = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    else:
        lambda_u, lambda_v, lambda_w = 3.0, 2.0, 1.0

    main(params_by_token = {'lambda_u': lambda_u, 'lambda_v': lambda_v, 'lambda_w': lambda_w})

