#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy

import matplotlib.pyplot as plt

from symodesys.firstorder import FirstOrderODESystem
from symodesys.integrator import SciPy_IVP_Integrator
#from symodesys.command_line import argument_parser


class VanDerPolOscillator(FirstOrderODESystem):

    num_dep_vars = 2
    num_params = 1

    # Following two lines are optional but useful for
    # automatic labeling when plotting:
    dep_var_tokens = 'u v'.split()
    # dep_var_tokens is used to generate sympy.Function()(indep_var) instances
    param_tokens = 'mu',
    #param_symbs   = sympy.symbols('mu,')

    @property
    def f(self):
        u, v = self.dep_var_func_symbs
        mu, = self.param_symbs
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
    main(params = {'mu': float(sys.argv[1])})
    # args = argument_parser.parse_args()
    # sys.exit(main(**vars(args)))
