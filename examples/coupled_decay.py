#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

import sympy
import numpy as np
import matplotlib.pyplot as plt

from symodesys.firstorder import SimpleFirstOrderODESystem
from symodesys.ivp import IVP

# TODO
# Implement automatic resolution of N - number of chained decays via
# Bateman's equations

class CoupledDecay(SimpleFirstOrderODESystem):

    # Following two lines are optional but useful for
    # automatic labeling when plotting:
    dep_var_tokens = 'u v w'.split()
    param_tokens   = 'lambda_u lambda_v lambda_w'.split()

    def init_f(self):
        u, v, w = self.dep_var_func_symbs
        lambda_u, lambda_v, lambda_w = self.param_symbs
        self.f = {u: -lambda_u * u,
                  v: lambda_u * u - lambda_v * v,
                  w: lambda_v * v - lambda_w * w,
                  }


    def analytic_u(self, indep_vals, y0):
        return y0[self['u']] * np.exp(-self.params_by_token['lambda_u']*indep_vals)


    def analytic_v(self, indep_vals, y0):
        return y0[self['v']] * np.exp(-self.params_by_token['lambda_v'] * indep_vals) + \
                 y0[self['u']] * self.params_by_token['lambda_u'] / \
                 (self.params_by_token['lambda_v'] - self.params_by_token['lambda_u']) * \
                 (np.exp(-self.params_by_token['lambda_u']*indep_vals) - \
                  np.exp( - self.params_by_token['lambda_v'] * indep_vals))

    def analytic_w(self, indep_vals, y0):
        return y0[self['w']] * np.exp(-self.params_by_token['lambda_w'] * indep_vals) + \
                 y0[self['v']] * self.params_by_token['lambda_v'] / \
                 (self.params_by_token['lambda_w'] - self.params_by_token['lambda_v']) * \
                 (np.exp(-self.params_by_token['lambda_v']*indep_vals) - \
                  np.exp(-self.params_by_token['lambda_w']*indep_vals)) + \
                 self.params_by_token['lambda_v'] * self.params_by_token['lambda_u'] * \
                 y0[self['u']] / (self.params_by_token['lambda_v'] - \
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

    t0 = 0.0
    tend = 1.5
    u0 = 7.0
    v0 = 5.0
    w0 = 3.0

    y0 = {cd['u']: u0,
          cd['v']: v0,
          cd['w']: w0,
          }

    ivp = IVP(cd, y0, t0)
    ivp._Integrator.abstol = 1e-8
    ivp._Integrator.reltol = 1e-8

    #N = 0 # adaptive stepsize controls output
    N = 100

    ivp.integrate(tend, N = N)

    t = ivp.tout
    analytic_u = cd.analytic_u(t, y0)
    analytic_v = cd.analytic_v(t, y0)
    analytic_w = cd.analytic_w(t, y0)

    uout = ivp.get_yout_by_symb(u)
    vout = ivp.get_yout_by_symb(v)
    wout = ivp.get_yout_by_symb(w)

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

    plt.plot(t, (uout - analytic_u) / ivp._Integrator.abstol,
             label = 'u abserr / abstol')
    plt.plot(t, (vout - analytic_v) / ivp._Integrator.abstol,
             label = 'v abserr / abstol')
    plt.plot(t, (wout - analytic_w) / ivp._Integrator.abstol,
             label = 'w abserr / abstol')
    plt.legend()

    plt.subplot(313)
    plt.plot(t, (uout - analytic_u) / analytic_u / ivp._Integrator.reltol,
             label = 'u relerr / reltol')
    plt.plot(t, (vout - analytic_v) / analytic_v / ivp._Integrator.reltol,
             label = 'v relerr / reltol')
    plt.plot(t, (wout - analytic_w) / analytic_w / ivp._Integrator.reltol,
             label = 'w relerr / reltol')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) > 1:
        lambda_u, lambda_v, lambda_w = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    else:
        lambda_u, lambda_v, lambda_w = 3.0, 2.0, 1.0

    main(params_by_token = {'lambda_u': lambda_u, 'lambda_v': lambda_v, 'lambda_w': lambda_w})

