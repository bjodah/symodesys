#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.firstorder import FirstOrderODESystem
from symodesys.integrator import SciPy_IVP_Integrator

from coupled_decay import CoupledDecay

# TODO
# Verify efficiency of general method compared to special use of
# Bateman's equations


def main(params_by_token):
    """
    """
    cd = CoupledDecay(params_by_token)

    u, v, w = cd.dep_var_func_symbs

    #N = 0 # adaptive stepsize controls output
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

    ivp = IVP(cd, y0)
    new_params = ivp.recurisive_analytic_reduction()

    intr = ivp.wrap_integrator(SciPy_IVP_Integrator)
    #intr = SciPy_IVP_Integrator(ivp._fo_odesys, new_params)

    int_kwargs = {'abstol': 1e-6,
                  'reltol': 1e-6}


    intr.integrate(y0, t0, tend, N, **int_kwargs)

    t = intr.tout

    uout = intr.get_yout_by_symb(u)
    vout = intr.get_yout_by_symb(v)
    wout = intr.get_yout_by_symb(w)

    plt.subplot(311)
    #intr.plot(interpolate = False, show = False)
    plt.plot(t, uout, '*', label = 'Numerical u')
    plt.plot(t, vout, 'o', label = 'Numerical v')
    plt.plot(t, wout, 'd', label = 'Numerical w')
    plt.plot(t, analytic_u, label = 'Analytic u')
    plt.plot(t, analytic_v, label = 'Analytic v')
    plt.plot(t, analytic_w, label = 'Analytic w')
    plt.legend()

    plt.subplot(312)

    plt.plot(t, (uout - analytic_u) / int_kwargs['abstol'],
             label = 'u abserr / abstol')
    plt.plot(t, (vout - analytic_v) / int_kwargs['abstol'],
             label = 'v abserr / abstol')
    plt.plot(t, (wout - analytic_w) / int_kwargs['abstol'],
             label = 'w abserr / abstol')
    plt.legend()

    plt.subplot(313)
    plt.plot(t, (uout - analytic_u) / analytic_u / int_kwargs['reltol'],
             label = 'u relerr / reltol')
    plt.plot(t, (vout - analytic_v) / analytic_v / int_kwargs['reltol'],
             label = 'v relerr / reltol')
    plt.plot(t, (wout - analytic_w) / analytic_w / int_kwargs['reltol'],
             label = 'w relerr / reltol')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        lambda_u, lambda_v, lambda_w = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    else:
        lambda_u, lambda_v, lambda_w = 3.0, 2.0, 1.0

    main(params_by_token = {'lambda_u': lambda_u, 'lambda_v': lambda_v, 'lambda_w': lambda_w})

