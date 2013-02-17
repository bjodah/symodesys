#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.firstorder import SimpleFirstOrderODESystem
from symodesys.integrator import SciPy_IVP_Integrator

# TODO
# Implement automatic resolution of N - number of chained decays via
# Bateman's equations

class CoupledDecay(SimpleFirstOrderODESystem):

    # Following two lines are optional but useful for
    # automatic labeling when plotting:
    dep_var_tokens = 'u v w'.split()
    param_tokens   = 'lambda_u lambda_v lambda_w'.split()

    def _init_f(self):
        u, v, w = self.dep_var_func_symbs
        lambda_u, lambda_v, lambda_w = self.param_symbs
        return {u: -lambda_u * u,
                v: lambda_u * u - lambda_v * v,
                w: lambda_v * v - lambda_w * w,
                }


    def analytic_u(self, indep_vals, y0):
        return y0['u'] * np.exp(-self.params_by_token['lambda_u']*indep_vals)


    def analytic_v(self, indep_vals, y0):
        return y0['v'] * np.exp(-self.params_by_token['lambda_v'] * indep_vals) + \
                 y0['u'] * self.params_by_token['lambda_u'] / \
                 (self.params_by_token['lambda_v'] - self.params_by_token['lambda_u']) * \
                 (np.exp(-self.params_by_token['lambda_u']*indep_vals) - \
                  np.exp( - self.params_by_token['lambda_v'] * indep_vals))

    def analytic_w(self, indep_vals, y0):
        return y0['w'] * np.exp(-self.params_by_token['lambda_w'] * indep_vals) + \
                 y0['v'] * self.params_by_token['lambda_v'] / \
                 (self.params_by_token['lambda_w'] - self.params_by_token['lambda_v']) * \
                 (np.exp(-self.params_by_token['lambda_v']*indep_vals) - \
                  np.exp(-self.params_by_token['lambda_w']*indep_vals)) + \
                 self.params_by_token['lambda_v'] * self.params_by_token['lambda_u'] * \
                 y0['u'] / (self.params_by_token['lambda_v'] - \
                            self.params_by_token['lambda_u']) * \
                 (1 / (self.params_by_token['lambda_w'] - \
                       self.params_by_token['lambda_u']) * \
                  (np.exp( - self.params_by_token['lambda_u'] * indep_vals) - \
                   np.exp( - self.params_by_token['lambda_w'] * indep_vals)) - \
                  1 / (self.params_by_token['lambda_w'] - \
                       self.params_by_token['lambda_v']) * \
                  (np.exp( - self.params_by_token['lambda_v'] * indep_vals) - \
                   np.exp( - self.params_by_token['lambda_w'] * indep_vals)))


def main(params_by_token):
    """
    """
    cd = CoupledDecay()
    cd.update_params_by_token(params_by_token)
    u, v, w = cd.dep_var_func_symbs
    intr = SciPy_IVP_Integrator(cd)
    int_kwargs = {'abstol': 1e-8,
                  'reltol': 1e-8}

    #N = 0 # adaptive stepsize controls output
    N = 100
    t0 = 0.0
    tend = 1.5
    u0 = 7.0
    v0 = 5.0
    w0 = 3.0

    y0 = {cd['u']: u0,
          cd['v']: v0,
          cd['w']: w0,
          }

    intr.integrate(y0, t0, tend, N, **int_kwargs)

    t = intr.tout
    analytic_u = u0 * np.exp(-lambda_u*t)
    analytic_v = v0 * np.exp(-lambda_v * t) + \
                 u0 * lambda_u / (lambda_v - lambda_u) * \
                 (np.exp(-lambda_u*t) - np.exp( - lambda_v * t))
    analytic_w = w0 * np.exp(-lambda_w * t) + \
                 v0 * lambda_v / (lambda_w - lambda_v) * \
                 (np.exp(-lambda_v*t) - np.exp(-lambda_w*t)) + \
                 lambda_v * lambda_u * u0 / (lambda_v - lambda_u) * \
                 (1 / (lambda_w - lambda_u) * (np.exp( - lambda_u * t) - \
                                               np.exp( - lambda_w * t)) - \
                  1 / (lambda_w - lambda_v) * (np.exp( - lambda_v * t) - \
                                               np.exp( - lambda_w * t)))

    uout = intr.get_yout_by_symb(u)
    vout = intr.get_yout_by_symb(v)
    wout = intr.get_yout_by_symb(w)

    plt.subplot(311)
    #intr.plot(interpolate = False, show = False)
    plt.plot(t, uout, '*', label = 'Numerical u')
    plt.plot(t, vout, 'o', label = 'Numerical v')
    plt.plot(t, wout, 'd', label = 'Numerical w')
    plt.plot(t, analytic_u, label = 'Analytic u')
    plt.plot(t, analytic_v, label = 'Analytic v')
    plt.plot(t, analytic_w, label = 'Analytic w')
    plt.legend()

    plt.subplot(312)

    plt.plot(t, (uout - analytic_u) / int_kwargs['abstol'],
             label = 'u abserr / abstol')
    plt.plot(t, (vout - analytic_v) / int_kwargs['abstol'],
             label = 'v abserr / abstol')
    plt.plot(t, (wout - analytic_w) / int_kwargs['abstol'],
             label = 'w abserr / abstol')
    plt.legend()

    plt.subplot(313)
    plt.plot(t, (uout - analytic_u) / analytic_u / int_kwargs['reltol'],
             label = 'u relerr / reltol')
    plt.plot(t, (vout - analytic_v) / analytic_v / int_kwargs['reltol'],
             label = 'v relerr / reltol')
    plt.plot(t, (wout - analytic_w) / analytic_w / int_kwargs['reltol'],
             label = 'w relerr / reltol')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        lambda_u, lambda_v, lambda_w = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    else:
        lambda_u, lambda_v, lambda_w = 3.0, 2.0, 1.0

    main(params_by_token = {'lambda_u': lambda_u, 'lambda_v': lambda_v, 'lambda_w': lambda_w})

