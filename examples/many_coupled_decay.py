#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.odesys import SimpleFirstOrderODESystem
from symodesys.ivp import IVP

"""
This example illustrates the ability of symodesys
to generate loop structures in code generation.
"""


class ManyCoupledDecay(SimpleFirstOrderODESystem):

    # Following two lines are optional but useful for
    # automatic labeling when plotting:
    dep_var_tokens = 'u v w'.split()
    param_tokens   = 'lambda_u lambda_v lambda_w'.split()

    @property
    def expressions(self):
        u, v, w = self['u'], self['v'], self['w']
        lambda_u, lambda_v, lambda_w = self.param_symbs
        return {u: -lambda_u * u,
                v: lambda_u * u - lambda_v * v,
                w: lambda_v * v - lambda_w * w,
                }

class DecayModel(DiffEqModel):

    # This should be able to relate to i-1 etc. ?
    param_tokens = 'lambda_${i}'

    @property
    def expression_model(self):
        return -self['lambda_${i}']*self.y


def main(params_by_token):
    """
    """
    mcd = ManyCoupledDecay()


if __name__ == '__main__':
    main()
