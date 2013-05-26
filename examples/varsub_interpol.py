#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP

from interpol import X5

# TODO: Let IVP do transfm and back-transfm internally + convert PiecewisePolynomial


def main(ODESys, y0, params, t0, tend, N = 0):
    """
    Integrate
    """

    # Setup
    odesys = ODESys()

    # Instantiate Initial Value Problem:
    ivp = IVP(odesys, y0, params, t0)

    # Set IVP to interally use variable substitution:
    # z=log(u)
    z = ivp.mk_depv_symbol('z')
    y = odesys['y']
    trnsfm = {z: sympy.log(y)}
    inv_trsfm = {y: sympy.exp(z)}

    # This call also converts initial values:
    ivp2 = ivp.use_internal_depv_trnsfm(trnsfm, inv_trsfm)

    # Solve IVP
    ivp2.integrate(tend, N=N)

    print '== t    ====='
    t = ivp2.integrator.tout
    print t
    print '== y(t) ===== (y=t**5/5+1.0, dydt = t**4))'
    num_y    = ivp2.trajectories()[ivp2.get_depv_from_token('y')][:,0]
    num_dydt = ivp2.trajectories()[ivp2.get_depv_from_token('y')][:,1]
    print num_y - (t**5/5+1.0)
    print num_dydt - t**4
    print '== z(t) ===== (z=log(y), dzdt = t**4*exp(-z))'
    num_z, num_dzdt = ivp2._Yres()[:,0,0], ivp2._Yres()[:,0,1]
    print num_z - np.log(t**5/5+1.0)
    print num_dzdt - t**4/(t**5/5+1.0)

    # Anlyse output
    plot_t = np.linspace(t0, tend, 50)

    ax = ivp2.plot(interpolate = True, show = False)
    for depv, cb in odesys.analytic_sol.items():
        ax.plot(plot_t, cb(odesys, plot_t, y0, params, t0), label='Analytic {}'.format(depv))

    plt.legend()
    plt.show()


if __name__ == '__main__':
    main(
        ODESys = X5,
        y0={'y': 1.0},
        t0=0.0,
        tend=3.0,
        params = {},
        N=10)
