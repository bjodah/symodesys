#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: Add interfaces to both GSL (CythonGSL or own wrapper?) and
# SUNDIALS (PySundials or own wrapper?)
# TODO: See if it is best to subclass sympy's codegen or use templates with Django template engine.

from __future__ import division

import warnings

import numpy as np
from scipy.interpolate import PiecewisePolynomial
import matplotlib.pyplot as plt

from symodesys.helpers import cache

class IVP_Integrator(object):
    """
    Integrates system of first order differential equations
    given initial values. Can be used directly with
    FirstOrderODESystem instance as fo_odesys, but often
    use of IVP_Problem is advised since it can do partial
    analytic treatment of problem at hand while preserving
    the notion of initial value even though some IV might
    have been reduced to parameters upon analytic reduction.
    """

    _dtype = np.float64

    def __init__(self, fo_odesys, params_by_symb = None):
        """

        Arguments:
        - `fo_odesys`: FO_ODESYS instance (initial value problem)
        - `params_by_symb`: Dictionary mapping parameter symbols to parameter values
        """
        self._fo_odesys = fo_odesys
        if params_by_symb == None: params_by_symb = {}
        for symb in params_by_symb.keys():
            if not symb in self._fo_odesys.param_symbs:
                raise KeyError('Parameter symbol {} not in ODE system'.format(symb))

        self._params_by_symb = {}
        for param_symb in self._fo_odesys.param_symbs:
            self._params_by_symb[param_symb] = params_by_symb.get(
                param_symb, self._fo_odesys.default_params[param_symb])
        self.post_init()

    def post_init(self):
        """
        Subclass for adding initialization steps to init
        """
        pass

    def update_params(self, params):
        self._params_by_symb.update(params)

    def compile(self):
        """
        To be subclassed.
        Should set self._compiled = True when done
        """
        pass

    def integrate(self, y0, t0, tend, N, abstol = None, reltol = None, h = None):
        """
        Should assign to self.tout and self.yout
        """
        pass

    def Dy(self):
        return np.array([self._fo_odesys.dydt(  t, self.yout[i,:], self.params_val_lst) for (i,), t \
                         in np.ndenumerate(self.tout)])

    def DDy(self):
        return np.array([self._fo_odesys.d2ydt2(t, self.yout[i,:], self.params_val_lst) for (i,), t \
                         in np.ndenumerate(self.tout)])

    @cache # never update tout, yout of an instance, create a new one instead
    def interpolators(self):
        Dy = self.Dy()
        DDy = self.DDy()
        intrpltrs = []
        for i in range(self.yout.shape[1]):
            intrpltrs.append(PiecewisePolynomial(
                self.tout, [[self.yout[j, i], Dy[j, i], DDy[j, i]] for j \
                            in range(self.yout.shape[0])], orders = 5))
        return intrpltrs

    def get_interpolated(self, t):
        return [self.interpolators()[i](t) for i in range(self.yout.shape[1])]

    def get_yout_by_symb(self, symb):
        return self.yout.view(
            dtype = [(str(x), self._dtype) for x \
                     in self._fo_odesys.dep_var_func_symbs])[str(symb)][:, 0]

    def plot(self, indices = None, interpolate = False, show = True):
        """
        Rudimentary plotting utility for quick inspection of solutions
        """
        if indices == None: indices = range(self.yout.shape[1])

        if interpolate:
            ipx = np.linspace(self.tout[0], self.tout[-1], 1000)
            ipy = np.array([self.get_interpolated(t) for t in ipx])
        ls = ['-', '--', ':']
        c = 'k b r g m'.split()
        m = 'o s t * d p h'.split()
        for i in indices:
            mi  = m[i % len(m)]
            lsi = ls[i % len(ls)]
            ci  = c[i % len(c)]
            lbl = str(self._fo_odesys.dep_var_func_symbs[i])
            if interpolate:
                plt.plot(ipx, ipy[:, i], label = lbl + ' (interpol.)',
                         marker = 'None', ls = lsi, color = ci)
                lsi = 'None'
            plt.plot(self.tout, self.yout[:, i], label = lbl,
                     marker = mi, ls = lsi, color = ci)
            plt.plot()
        plt.legend()
        if show: plt.show()

    def init_yout_tout_for_fixed_step_size(self, t0, tend, N):
        dt = (tend - t0) / (N - 1)
        self.tout = np.linspace(t0, tend, N)
        # Handle other dtype for tout here? linspace doesn't support dtype arg..
        self.yout = np.zeros((N, self._fo_odesys.num_dep_vars), dtype = self._dtype)

    @property
    def params_val_lst(self):
        return [self._params_by_symb[k] for k in self._fo_odesys.param_symbs]

class SciPy_IVP_Integrator(IVP_Integrator):

    def post_init(self):
        from scipy.integrate import ode
        self._r = ode(self._fo_odesys.dydt, self._fo_odesys.dydt_jac)


    def integrate(self, y0, t0, tend, N, abstol = None, reltol = None, h = None):
        y0_val_lst = [y0[k] for k in self._fo_odesys.dep_var_func_symbs]
        self._r.set_initial_value(y0_val_lst)
        assert len(self._params_by_symb) == self._fo_odesys.num_params
        self._r.set_f_params(self.params_val_lst)
        self._r.set_jac_params(self.params_val_lst)
        if N > 0:
            # Fixed stepsize
            self._r.set_integrator('vode', method = 'bdf', with_jacobian = True)
            self.init_yout_tout_for_fixed_step_size(t0, tend, N)
            for i, t in enumerate(self.tout):
                if i == 0:
                    self.yout[i, :] = y0_val_lst
                    continue
                self.yout[i, :] = self._r.integrate(self.tout[i])
                assert self._r.successful()
        else:
            # Adaptive step size reporting
            # http://stackoverflow.com/questions/12926393/\
            #   using-adaptive-step-sizes-with-scipy-integrate-ode
            self._r.set_integrator('vode', method = 'bdf', with_jacobian = True, nsteps = 1)
            self._r._integrator.iwork[2] =- 1
            yout, tout= [], []
            warnings.filterwarnings("ignore", category=UserWarning)
            yout.append(y0_val_lst)
            tout.append(t0)
            while self._r.t < tend:
                self._r.integrate(tend, step=True)
                yout.append(self._r.y)
                tout.append(self._r.t)
            warnings.resetwarnings()
            self.yout = np.array(yout)
            self.tout = np.array(tout)


class Mpmath_IVP_Integrator(IVP_Integrator):
    """
    Only for demonstration purposes
    """

    def post_init(self):
        pass

    def integrate(self, y0, t0, tend, N, abstol = None, reltol = None, h = None):
        y0_val_lst = [y0[k] for k in self._fo_odesys.dep_var_func_symbs]
        cb = lambda x, y: self._fo_odesys.dydt(x, y, self.params_val_lst)
        self._num_y = sympy.mpmath.odefun(cb, t0, y0_val_lst, tol = abstol)
        if N > 0:
            # Fixed stepsize
            self.init_yout_tout_for_fixed_step_size(t0, tend, N)
            for i, t in enumerate(self.tout):
                if i == 0:
                    self.yout[i, :] = y0_val_lst
                    continue
                self.yout[i, :] = self._num_y(self.tout[i])
        else:
            self.yout = np.array([y0_val_lst, self._num_y(tend)])
            self.tout = np.array([t0, tend])


class GSL_IVP_Integrator(IVP_Integrator):
    pass
