#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from decay import Decay


def main(ODESys, y0, t0, tend, params, N=0):
    """
    This example does not use the transform hiding mechanisms
    of IVP but directly acts on FirstOrderODESystem methods
    for variable transformations.

    For a more convinient way  to do variable substitutions look
    at varsub_interpol.py
    """
    odesys = ODESys()

    # Substitute u for v=log(u)
    z = sympy.Function('z')(odesys.indepv)
    u = odesys['u']
    trnsfm = {z: sympy.log(u)}
    inv_trsfm = {u: sympy.exp(z)}
    new_odesys = odesys.transform_depv(trnsfm, inv_trsfm)
    sympy.pprint(odesys.eqs)
    # Convert initial values:
    new_y0 = {new_odesys[dv]: y0[dv] if dv in y0 else trnsfm[dv].subs(
        odesys.ensure_dictkeys_as_symbs(y0)) for \
          dv in new_odesys.all_depv}

    ivp = IVP(new_odesys, new_y0, params, t0)
    ivp.integrate(tend, N = N)
    t, z = ivp.indep_out(), ivp.trajectories()[new_odesys['z']]

    ax = plt.subplot(311)
    ivp.plot(ax=ax, interpolate=True, show=False)
    ax.legend()
    # ax.xlabel('t')
    # ax.ylabel('z')

    analytic_u = odesys.analytic_u(t, y0, params, t0)
    analytic_z = np.log(analytic_u)

    plt.subplot(312)
    plt.plot(t, z[:, 0] - analytic_z,
             label = 'abserr')
    plt.legend()
    plt.subplot(313)
    plt.plot(t, (z[:, 0] - analytic_z) / analytic_z,
             label = 'relerr')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    main(
        Decay,
        y0 = {'u': 1.0},
        t0 = 0.0,
        tend = 10.0,
        params = {'lambda_u': 0.2},
        N = 10)
