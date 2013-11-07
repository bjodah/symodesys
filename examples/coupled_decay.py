#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np

from symodesys.odesys import SimpleFirstOrderODESystem

# TODO
# Implement automatic resolution of N - number of chained decays via
# Bateman's equations (Note: breaks down when lambda_k = lambda_l)

class CoupledDecay(SimpleFirstOrderODESystem):
    """
    A first order system of 3 ODEs (u, v, w)
    with 3 parameters (lambda_u, lambda_v, lambda_w)
    `analytic_sol` contains the analytic solutions (valid for
    lambda_u != lambda_v != lambda_w)
    """
    # Following two lines are optional but useful for
    # automatic labeling when plotting:
    depv_tokens = 'u v w'.split()
    param_tokens   = 'lambda_u lambda_v lambda_w'.split()

    @property
    def expressions(self):
        u, v, w = self['u'], self['v'], self['w']
        lambda_u, lambda_v, lambda_w = self.param_symbs
        return {u: -lambda_u * u,
                v: lambda_u * u - lambda_v * v,
                w: lambda_v * v - lambda_w * w,
                }

    def analytic_u(self, indep_vals, y0, params, t0):
        return y0['u'] * np.exp(-params['lambda_u']*indep_vals)


    def analytic_v(self, indep_vals, y0, params, t0):
        return y0['v'] * np.exp(-params['lambda_v'] * indep_vals) + \
                 y0['u'] * params['lambda_u'] / \
                 (params['lambda_v'] - params['lambda_u']) * \
                 (np.exp(-params['lambda_u']*indep_vals) - \
                  np.exp( - params['lambda_v'] * indep_vals))

    def analytic_w(self, indep_vals, y0, params, t0):
        return y0['w'] * np.exp(-params['lambda_w'] * indep_vals) + \
                 y0['v'] * params['lambda_v'] / \
                 (params['lambda_w'] - params['lambda_v']) * \
                 (np.exp(-params['lambda_v']*indep_vals) - \
                  np.exp(-params['lambda_w']*indep_vals)) + \
                 params['lambda_v'] * params['lambda_u'] * \
                 y0['u'] / (params['lambda_v'] - \
                            params['lambda_u']) * \
                 (1 / (params['lambda_w'] - \
                       params['lambda_u']) * \
                  (np.exp( - params['lambda_u'] * indep_vals) - \
                   np.exp( - params['lambda_w'] * indep_vals)) - \
                  1 / (params['lambda_w'] - \
                       params['lambda_v']) * \
                  (np.exp( - params['lambda_v'] * indep_vals) - \
                   np.exp( - params['lambda_w'] * indep_vals)))

    analytic_sol = {'u': analytic_u, 'v': analytic_v, 'w': analytic_w}


def coupled_decay_numeric_vs_analytic(
        Integrator=None, integrator_kwargs=None, logger=None,
        **kwargs):
    if Integrator == None:
        from symodesys.integrator import SciPy_IVP_Integrator
        Integrator = SciPy_IVP_Integrator
    integrator_kwargs = integrator_kwargs or {}

    if not 'acceptance_factor' in kwargs: kwargs['acceptance_factor'] = 100
    from symodesys import numeric_vs_analytic
    numeric_vs_analytic(
        ODESys = CoupledDecay,
        depv_init = {'u': 1.0, 'v': 1.0, 'w': 1.0},
        params = {'lambda_u': 1/3,'lambda_v': 1/5,'lambda_w': 1/7},
        indepv_init = 0,
        indepv_end = 5.0,
        integrator = Integrator(**integrator_kwargs),
        N = 100,
        plot=True,
        logger=logger,
        **kwargs
    )


if __name__ == '__main__':
    coupled_decay_numeric_vs_analytic()
