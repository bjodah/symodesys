#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

from symodesys.odesys import SimpleFirstOrderODESystem
from symodesys.ivp import IVP


class VanDerPolOscillator(SimpleFirstOrderODESystem):

    dep_var_tokens = 'u v'.split()
    param_tokens = 'mu',

    @property
    def expressions(self):
        u, v, mu = [self[x] for x in 'u v mu'.split()]
        return {u: v,
                v: -u + mu*v*(1 - u**2),
                }


def main(y0, params, tend, t0 = 0.0, N = 50):
    """
    Example program integrating an IVP problem of van der Pol oscillator
    default is adaptive step size (N=0)
    """
    vdpo = VanDerPolOscillator()
    ivp = IVP(vdpo, y0, params, t0)

    ivp.integrate(tend, N)
    ivp.plot(interpolate = True, datapoints=False, show = True)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        mu = float(sys.argv[1])
    else:
        mu = 1.0

    main(y0={'u': 1.0, 'v': 0.0}, params = {'mu': mu}, tend=10.0)
