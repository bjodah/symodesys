#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import logging

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.ivp import IVP
from symodesys.integrator import Mpmath_IVP_Integrator
from decay import Decay


def main(Sys, y0, t0, tend, params, N=20, nderiv=2, logger=None):
    """
    Integrate
    """
    odesys = Sys()

    ivp = IVP(odesys, y0, params, t0,
              integrator=Mpmath_IVP_Integrator(nderiv=nderiv),
              logger=logger)

    # Attempt analytic reduction
    ivp.recursive_analytic_reduction(complexity=3)
    ivp.integrate(tend, N=N)
    t, u = ivp.indepv_out(), ivp.trajectories()[odesys['u']]

    expected_u = np.vstack((
        y0['u']*np.exp(-params['lambda_u']*t),
        y0['u']*np.exp(-params['lambda_u']*t)*-params['lambda_u'],
        y0['u']*np.exp(-params['lambda_u']*t)*params['lambda_u']**2,
        )).transpose()
    assert np.allclose(u, expected_u[:,:nderiv+1])

    # Manuallh make a plot
    ax=plt.subplot(311)
    ivp.plot(ax=ax)

    analytic_u = odesys.analytic_u(t, y0, params, t0)
    plt.subplot(312)
    plt.plot(t, u[:, 0] - analytic_u, label = 'abserr')
    plt.legend()
    plt.subplot(313)
    plt.plot(t, (u[:, 0] - analytic_u) / analytic_u, label = 'relerr')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    logger = logging.getLogger(__file__)
    main(
        Sys = Decay,
        y0={'u': 1.0},
        t0=0.0,
        tend=3.0,
        params={'lambda_u': 0.2},
        logger=logger
        )
