#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from decay import Decay


def main_analytic(ODESys, y0, t0, tend, params N = 0):
    """
    This example does not use the transform hiding mechanisms
    of IVP but directly acts on FirstOrderODESystem methods
    for variable transformations.
    """
    odesys = ODESys()

    # Substitute u for v=log(u)
    z = sympy.Function('z')(odesys.indepv)
    u = odesys['u']
    trnsfm = {z: sympy.log(u)}
    inv_trsfm = {u: sympy.exp(z)}
    odesys = odesys.transform_depv(trnsfm, inv_trsfm)
    sympy.pprint(odesys.eqs)
    # Convert initial values:
    y0 = {odesys[dv]: y0[dv] if dv in y0 else trnsfm[dv].subs(y0) for \
          dv in odesys.all_depv}

    ivp = IVP(odesys, y0, param_vals_by_symb, t0)
    ivp.integrate(tend, N = N)
    t, z = ivp.tout, ivp.yout

    plt.subplot(311)
    ivp.plot(interpolate = True, show = False)

    analytic_u = odesys.analytic_u(
        t, init_dep_var_vals_by_token, param_vals_by_symb)
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
    plot_numeric_vs_analytic(
        Decay,
        y0 = {'u': 1.0},
        t0 = 0.0,
        tend = 10.0,
        params = {'lambda_u': 0.2},
        N = 0)
