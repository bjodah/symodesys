#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy

from symodesys.firstorder import FirstOrderODESystem
from symodesys.integrator import SciPy_IVP_Integrator

class VanDerPolOscillator(FirstOrderODESystem):

    dep_var_tokens = 'u v'.split()
    param_tokens = 'mu',

    @property
    def f(self):
        u, v, mu = [self[x] for x in 'u v mu'.split()]
        return {u: v,
                v: -u + mu*v*(1 - u**2),
                }


def main(params):
    """
    """
    vdpo = VanDerPolOscillator()
    vdpo_params = dict([(vdpo[k], v) for k, v in params.iteritems()])
    intr = SciPy_IVP_Integrator(vdpo, vdpo_params)

    y0 = {vdpo['u']: 1.0, vdpo['v']: 0.0}

    int_kwargs = {'abstol': 1e-6,
                  'reltol': 1e-6}

    #N = 0 # adaptive stepsize controls output
    N = 100
    intr.integrate(y0, 0.0, 10.0, N, **int_kwargs)
    print intr.tout.shape, intr.yout.shape
    intr.plot(interpolate = True)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        mu = float(sys.argv[1])
    else:
        mu = 1.0

    main(params = {'mu': mu})

    # args = argument_parser.parse_args()
    # sys.exit(main(**vars(args)))
