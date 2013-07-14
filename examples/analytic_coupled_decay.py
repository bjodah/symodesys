#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging

import matplotlib.pyplot as plt

from symodesys.ivp import IVP

from coupled_decay import CoupledDecay


def main(params, logger=None):
    """
    Solve a couled decay network analytically
    using sympy's symbolic algortihms and compare
    with hand derived solution.
    """
    cd = CoupledDecay()

    u, v, w = cd['u'], cd['v'], cd['w']

    N = 100
    t0 = 0.0
    tend = 1.5

    y0 = {'u': 7,
          'v': 5,
          'w': 3,
          }

    ivp = IVP(cd, y0, params, t0, logger=logger)
    new_params = ivp.recursive_analytic_reduction(complexity=3)

    ivp.integrator.abstol = 1e-6
    ivp.integrator.reltol = 1e-6

    ivp.integrate(tend, N = N)

    t = ivp.indepv_out()
    analytic_u = cd.analytic_u(t, y0, params, t0)
    analytic_v = cd.analytic_v(t, y0, params, t0)
    analytic_w = cd.analytic_w(t, y0, params, t0)

    traj = ivp.trajectories()
    uout, vout, wout = [traj[cd[x]][:,0] for x in 'u v w'.split()]

    plt.subplot(311)
    plt.plot(t, uout, '*', label = 'Numerical u')
    plt.plot(t, vout, 'o', label = 'Numerical v')
    plt.plot(t, wout, 'd', label = 'Numerical w')
    plt.plot(t, analytic_u, label = 'Analytic u')
    plt.plot(t, analytic_v, label = 'Analytic v')
    plt.plot(t, analytic_w, label = 'Analytic w')
    plt.legend()

    plt.subplot(312)
    plt.plot(t, (uout - analytic_u) / ivp.integrator.abstol,
             label = 'u abserr / abstol')
    plt.plot(t, (vout - analytic_v) / ivp.integrator.abstol,
             label = 'v abserr / abstol')
    plt.plot(t, (wout - analytic_w) / ivp.integrator.abstol,
             label = 'w abserr / abstol')
    plt.legend()

    plt.subplot(313)
    plt.plot(t, (uout - analytic_u) / analytic_u / ivp.integrator.reltol,
             label = 'u relerr / reltol')
    plt.plot(t, (vout - analytic_v) / analytic_v / ivp.integrator.reltol,
             label = 'v relerr / reltol')
    plt.plot(t, (wout - analytic_w) / analytic_w / ivp.integrator.reltol,
             label = 'w relerr / reltol')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        lambda_u, lambda_v, lambda_w = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    else:
        lambda_u, lambda_v, lambda_w = 3.0, 2.0, 1.0

    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__file__)
    main(params={'lambda_u': lambda_u, 'lambda_v': lambda_v, 'lambda_w': lambda_w},
         logger=logger)
