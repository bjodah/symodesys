#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.odesys import SimpleFirstOrderODESystem
from symodesys.ivp import IVP
from symodesys.gsl import GSL_IVP_Integrator

"""
This example illustrates the need of symodesys
to generate loop structures in code generation.
(and efficiently handle sparse systems)
"""


class BirthDeathSystem(SimpleFirstOrderODESystem):

    # Following two lines are optional but useful for
    # automatic labeling when plotting:

    @classmethod
    def mk_of_dim(cls, dim):
        """
        Makes an instance of provided `dim`
        """
        cls.dep_var_tokens = ['y{}'.format(i) for i in range(dim)]
        cls.param_tokens   = ['l{}'.format(i) for i in range(dim)]
        cls.dim = dim
        return cls()

    @property
    def expressions(self):
        y = [self['y{}'.format(i)] for i in range(self.dim)]
        l = [self['l{}'.format(i)] for i in range(self.dim)]
        first_last = {y[0]: -l[0]*y[0],
                      y[self.dim-1]: -l[self.dim-2]*y[self.dim-1]}
        first_last.update({
            y[i]: y[i-1]*l[i-1]-y[i]*l[i] for \
            i in range(1,self.dim-1)
        })
        return first_last

def main(n):
    """
    Plot the evolutoin of a birth death network
    """

    bd = BirthDeathSystem.mk_of_dim(n)
    y0 = {'y{}'.format(i): (n-i)*1.0 for i in range(n)}
    params = {'l{}'.format(i): (n-i)/10.0 for i in range(n)}
    t0, tend, N = 0, 0.5, 100
    ivp = IVP(bd, y0, params, t0, integrator=GSL_IVP_Integrator())
    trnsfm, inv_trnsfm = {}, {}
    for i in range(n):
        # Generate the transform (incl. inverse)
        lny = sympy.Symbol('lny{}'.format(i))
        y = bd['y{}'.format(i)]
        trnsfm[lny] = sympy.log(y)
        inv_trnsfm[y] = sympy.exp(lny)
    # Apply the transform
    ivp.use_internal_depv_trnsfm(trnsfm, inv_trnsfm)
    ivp.integrate(tend, N)
    ivp.plot(interpolate=True, datapoints=False, show=True)


if __name__ == '__main__':
    main(10)
