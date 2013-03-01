#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy

from symodesys.firstorder import SimpleFirstOrderODESystem
from symodesys.ivp import IVP


class VanDerPolOscillator(SimpleFirstOrderODESystem):

    dep_var_tokens = 'u v'.split()
    param_tokens = 'mu',

    def init_f(self):
        u, v, mu = [self[x] for x in 'u v mu'.split()]
        self.f = {u: v,
                  v: -u + mu*v*(1 - u**2),
                  }


def main(params):
    """
    Example program integrating an IVP problem of van der Pol oscillator
    """
    vdpo = VanDerPolOscillator()
    vdpo.update_params_by_token(params)

    y0 = {vdpo['u']: 1.0, vdpo['v']: 0.0}
    t0 = 0.0
    ivp = IVP(vdpo, y0, t0)

    # TODO: add abstraction layer for _Integrator.abstol etc?
    ivp._Integrator.abstol = 1e-6
    ivp._Integrator.reltol = 1e-6

    N = 100
    ivp.integrate(10.0, N)
    ivp.plot(interpolate = True, show = True)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        mu = float(sys.argv[1])
    else:
        mu = 1.0

    main(params = {'mu': mu})

