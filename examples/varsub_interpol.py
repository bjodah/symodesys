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
    z = ivp.fo_odesys.mk_depv('z')
    y = odesys['y']
    trnsfm = {z: sympy.log(y)}
    inv_trsfm = {y: sympy.exp(z)}

    # This call also converts initial values:
    ivp2 = ivp.use_internal_depv_trnsfm(trnsfm, inv_trsfm)

    # Solve IVP
    ivp2.integrate(tend, N=N)

    num_y    = ivp2.trajectories()[ivp2.get_depv_from_token('y')][:,0]
    num_dydt = ivp2.trajectories()[ivp2.get_depv_from_token('y')][:,1]
    t = ivp2.indepv_out()
    assert np.allclose(num_y, (t**5/5+1.0))
    assert np.allclose(num_dydt, t**4)

    num_z, num_dzdt = ivp2._Yres()[:,0,0], ivp2._Yres()[:,0,1]
    assert np.allclose(num_z, np.log(t**5/5+1.0))
    assert np.allclose(num_dzdt, t**4/(t**5/5+1.0))

    # Anlyse output
    plot_t = np.linspace(t0, tend, 50)

    ax = ivp2.plot(interpolate = True, show = False, ax=plt.subplot(2, 1, 1))
    for depv, cb in odesys.analytic_sol.items():
        ax.plot(plot_t, cb(odesys, plot_t, y0, params, t0), label='Analytic {}'.format(depv))
    plt.legend()

    plt.subplot(2, 1, 2)

    for depv, cb in odesys.analytic_sol.items():
        analytic = cb(odesys, ivp2.indepv_out(),
                      y0, params, t0)
        numeric = ivp2.trajectories()[odesys[depv]][:,0]
        plt.plot(ivp2.indepv_out(), analytic-numeric, label='Error {}'.format(depv))
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
