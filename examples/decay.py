#!/usr/bin/env python
# -*- coding: utf-8 -*-

# External imports
import numpy as np

# Package imports
from symodesys import SimpleFirstOrderODESystem
from symodesys import numeric_vs_analytic

def analytic_decay(indep_vals, y0, params, t0):
    return y0['u'] * np.exp(-params['lambda_u']*indep_vals)

class Decay(SimpleFirstOrderODESystem):

    depv_tokens = 'u',
    param_tokens   = 'lambda_u',

    @property
    def expressions(self):
        return {self['u']: self['lambda_u'] * -self['u']}

    @property
    def analytic_sol(self):
        return {self['u']: analytic_decay}


if __name__ == '__main__':
    numeric_vs_analytic(
        Decay, {'u': 1.0}, {'lambda_u': 0.2}, 0.0,
        10.0, N=30, acceptance_factor=1e3)
