#!/usr/bin/env python
# -*- coding: utf-8 -*-

# External imports
import numpy as np

# Package imports
from symodesys import SimpleFirstOrderODESystem
from symodesys.convenience import plot_numeric_error

class Decay(SimpleFirstOrderODESystem):

    depv_tokens = 'u',
    param_tokens   = 'lambda_u',

    @property
    def expressions(self):
        return {self['u']: self['lambda_u'] * -self['u']}

    def analytic_u(self, indep_vals, y0, params, t0):
        return y0['u'] * np.exp(-params['lambda_u']*indep_vals)

    analytic_sol = {'u': analytic_u}


if __name__ == '__main__':
    plot_numeric_error(Decay, {'u': 1.0}, {'lambda_u': 0.2}, 0.0, 10.0, N=30)
