#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: Add interfaces to both GSL (CythonGSL or own wrapper?) and
# SUNDIALS (PySundials or own wrapper?)
# TODO: See if it is best to subclass sympy's codegen or use templates with Django template engine.

from __future__ import division

import warnings

import sympy.mpmath
import numpy as np


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

    abstol = 1e-6 # Default absolute tolerance
    reltol = 1e-6 # Default relative tolerance

    kwargs_attr=['abstol', 'reltol']

    def __init__(self, **kwargs):
        """

        Arguments:
        - `fo_odesys`: FO_ODESYS instance (initial value problem)
        """
        self._fo_odesys = None
        for attr in kwargs:
            if attr in self.kwargs_attr:
                setattr(self, attr, kwargs.pop(attr))
            else:
                raise AttributeError('Unkown kwarg: {}'.format(attr))

    def set_fo_odesys(self, fo_odesys):
        self._fo_odesys = fo_odesys

    def run(self, y0, t0, tend, param_vals,
                  N, abstol = None, reltol = None, h = None,
                  order = 0):
        """
        Should assign to self.tout and self.yout
        - `y0`: Dict mapping dep. var. symbs to initial values
        - `t0`: Floating point value of initial value for indep var
        - `tend`: Integration limit for the IVP (max val of indep var)
        - `N`: number of output values. N = 0 signals adaptive step size and dynamic
               allocation of length
        - `abstol`: absolute tolerance for numerical integrator driver routine
        - `reltol`: relative tolerance for numerical integrator driver routine
        - `order`: Up to what order should d^n/dt^n derivatives of y be evaluated
                   (defaults to 0)
        Changes to the signature of this function must be propagated to IVP.integrate
        """
        pass


    def init_Yout_tout_for_fixed_step_size(self, t0, tend, N, order = 0):
        dt = (tend - t0) / (N - 1)
        NY = len(self._fo_odesys.non_analytic_depv)
        self.tout = np.asarray(np.linspace(t0, tend, N), dtype = self._dtype)
        # Handle other dtype for tout here? linspace doesn't support dtype arg..
        self.Yout = np.zeros((N, NY, order + 1), dtype = self._dtype)


class Mpmath_IVP_Integrator(IVP_Integrator):
    """
    Only for demonstration purposes
    """

    def run(self, y0, t0, tend, param_vals,
                  N, abstol = None, reltol = None, h = None,
            order=0):
        y0_val_lst = [y0[k] for k in self._fo_odesys.non_analytic_depv]
        param_val_lst = self._fo_odesys.param_val_lst(param_vals)
        cb = lambda x, y: self._fo_odesys.dydt(x, y, param_val_lst)
        self._num_y = sympy.mpmath.odefun(cb, t0, y0_val_lst, tol = abstol)
        if N > 0:
            # Fixed stepsize
            self.init_Yout_tout_for_fixed_step_size(t0, tend, N, order)
            for i, t in enumerate(self.tout):
                if i == 0:
                    self.Yout[i, :, 0] = y0_val_lst
                else:
                    self.Yout[i, :, 0] = self._num_y(self.tout[i])

                if order > 0:
                    self.Yout[i, :, 1] = self._fo_odesys.dydt(
                        self.tout[i], self.Yout[i, :, 0], param_val_lst)
                if order > 1:
                    self.Yout[i, :, 2] = self._fo_odesys.d2ydt2(
                        self.tout[i], self.Yout[i, :, 0], param_val_lst)

        else:
            raise NotImplementedError('Mpmath_IVP_Integrator does not currently support'+\
                                      ' adaptive step-size control, please consider using'+\
                                      ' SciPy_IVP_Integrator if that is what is sought for.')

            # NOTE: It is not quite clear to me how to extract
            # intermediate steps taken by mpmath.odefun... hence NotImplementedError
            # self.tout = np.array([t0, tend])
            # # Set Yout depending on order
            # yout = np.array([y0_val_lst, self._num_y(tend)])
            # if order == 0:
            #     self.Yout = yout.reshape(yout.shape+(1,))
            # if order > 0:
            #     self.Yout = np.concatenate((yout, np.array(
            #         [self._fo_odesys.dydt(
            #             self.tout[i], yout[i, :], param_val_lst) for i in range(len(self.tout))])), axis=2)
            # if order > 1:
            #     self.Yout = np.concatenate((yout, np.array(
            #         [self._fo_odesys.d2ydt2(
            #             self.tout[i], yout[i, :], param_val_lst) for i in range(len(self.tout))])), axis=2)
            # self.Yout = yout


class SciPy_IVP_Integrator(IVP_Integrator):
    """
    Uses integrator provided in SciPy,
    beware that this method is slow (expensive sympy.subs calls during
    evaluation of jacobian and dydt), consider using a native code
    generating integrator such as e.g. GSL_IVP_Integrator
    """

    def __init__(self, **kwargs):
        super(SciPy_IVP_Integrator, self).__init__(**kwargs)

    def set_fo_odesys(self, fo_odesys):
        super(SciPy_IVP_Integrator, self).set_fo_odesys(fo_odesys)
        from scipy.integrate import ode
        self._r = ode(self._fo_odesys.dydt, self._fo_odesys.dydt_jac)

    def run(self, y0, t0, tend, param_vals,
                  N, abstol=None, reltol=None, h=None,
                  order=0):
        y0_val_lst = [y0[k] for k in self._fo_odesys.non_analytic_depv]
        param_val_lst = self._fo_odesys.param_val_lst(param_vals)
        self._r.set_initial_value(y0_val_lst, t0)
        self._r.set_f_params(param_val_lst)
        self._r.set_jac_params(param_val_lst)
        if N > 0:
            # Fixed stepsize
            self._r.set_integrator('vode', method = 'bdf', with_jacobian = True)
            self.init_Yout_tout_for_fixed_step_size(t0, tend, N, order)
            for i, t in enumerate(self.tout):
                if i == 0:
                    self.Yout[i, :, 0] = y0_val_lst
                else:
                    self.Yout[i, :, 0] = self._r.integrate(self.tout[i])

                if order > 0:
                    self.Yout[i, :, 1] = self._fo_odesys.dydt(
                        self.tout[i], self.Yout[i, :, 0], param_val_lst)
                if order > 1:
                    self.Yout[i, :, 2] = self._fo_odesys.d2ydt2(
                        self.tout[i], self.Yout[i, :, 0], param_val_lst)
                assert self._r.successful()
        else:
            # Adaptive step size reporting
            # http://stackoverflow.com/questions/12926393/\
            #   using-adaptive-step-sizes-with-scipy-integrate-ode
            self._r.set_integrator('vode', method = 'bdf', with_jacobian = True, nsteps = 1)
            self._r._integrator.iwork[2] =- 1
            tout, yout, dyout, ddyout = [], [], [], []
            warnings.filterwarnings("ignore", category=UserWarning)
            # yout.append(np.array(y0_val_lst))
            # tout.append(t0)
            keep_going = True
            while keep_going:
                keep_going = self._r.t < tend
                if order > 0:
                    dyout.append(self._fo_odesys.dydt(
                        self._r.t, self._r.y, param_val_lst))
                if order > 1:
                    ddyout.append(self._fo_odesys.d2ydt2(
                        self._r.t, self._r.y, param_val_lst))
                self._r.integrate(tend, step=True)
                yout.append(self._r.y)
                tout.append(self._r.t)
            warnings.resetwarnings()
            tout = np.array(tout)
            yout, dyout, ddyout = map(np.array, (yout, dyout, ddyout))
            if order == 0:
                self.Yout = yout.reshape((yout.shape[0], yout.shape[1], 1))
            elif order == 1:
                self.Yout = np.concatenate((yout[...,np.newaxis], dyout[...,np.newaxis]), axis=2)
            else:
                self.Yout = np.concatenate((yout[...,np.newaxis], dyout[...,np.newaxis],
                                           ddyout[...,np.newaxis]), axis=2)
            self.tout = tout
