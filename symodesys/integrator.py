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
    nderiv = 1 # Default maximum order of derivatives

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
                  N, abstol=None, reltol=None, h=None,
                  nderiv=None):
        """
        Should assign to self.tout and self.yout
        - `y0`: Dict mapping dep. var. symbs to initial values
        - `t0`: Floating point value of initial value for indep var
        - `tend`: Integration limit for the IVP (max val of indep var)
        - `N`: number of output values. N = 0 signals adaptive step size and dynamic
               allocation of length
        - `abstol`: absolute tolerance for numerical integrator driver routine
        - `reltol`: relative tolerance for numerical integrator driver routine
        - `nderiv`: Up to what nderiv should d^n/dt^n derivatives of y be evaluated
                   (defaults to 0)
        Changes to the signature of this function must be propagated to IVP.integrate
        """
        pass

    def clean(self):
        """
        Clean up temporary files.
        """
        pass

    def init_Yout_tout_for_fixed_step_size(self, t0, tend, N):
        dt = (tend-t0) / (N-1)
        NY = len(self._fo_odesys.non_analytic_depv)
        self.tout = np.asarray(np.linspace(t0, tend, N), dtype = self._dtype)
        # Handle other dtype for tout here? linspace doesn't support dtype arg..
        self.Yout = np.zeros((N, NY, self.nderiv+1), dtype = self._dtype)


class Mpmath_IVP_Integrator(IVP_Integrator):
    """
    Only for demonstration purposes
    """

    def run(self, y0, t0, tend, param_vals,
                  N, abstol = None, reltol = None, h = None,
            nderiv=None):
        self.nderiv = nderiv or self.nderiv
        y0_val_lst = [y0[k] for k in self._fo_odesys.non_analytic_depv]
        param_val_lst = self._fo_odesys.param_val_lst(param_vals)
        cb = lambda x, y: self._fo_odesys.dydt(x, y, param_val_lst)
        self._num_y = sympy.mpmath.odefun(cb, t0, y0_val_lst, tol = abstol)
        if N > 0:
            # Fixed stepsize
            self.init_Yout_tout_for_fixed_step_size(t0, tend, N)
            for i, t in enumerate(self.tout):
                if i == 0:
                    self.Yout[i, :, 0] = y0_val_lst
                else:
                    self.Yout[i, :, 0] = self._num_y(self.tout[i])

                if self.nderiv > 0:
                    self.Yout[i, :, 1] = self._fo_odesys.dydt(
                        self.tout[i], self.Yout[i, :, 0], param_val_lst)
                if self.nderiv > 1:
                    self.Yout[i, :, 2] = self._fo_odesys.d2ydt2(
                        self.tout[i], self.Yout[i, :, 0], param_val_lst)

        else:
            raise NotImplementedError('Mpmath_IVP_Integrator does not currently support'+\
                                      ' adaptive step-size control, please consider using'+\
                                      ' SciPy_IVP_Integrator if that is what is sought for.')


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
                  nderiv=None):
        self.nderiv = nderiv or self.nderiv
        y0_val_lst = [y0[k] for k in self._fo_odesys.non_analytic_depv]
        param_val_lst = self._fo_odesys.param_val_lst(param_vals)
        self._r.set_initial_value(y0_val_lst, t0)
        self._r.set_f_params(param_val_lst)
        self._r.set_jac_params(param_val_lst)
        if N > 0:
            # Fixed stepsize
            self._r.set_integrator('vode', method = 'bdf', with_jacobian = True)
            self.init_Yout_tout_for_fixed_step_size(t0, tend, N)
            for i, t in enumerate(self.tout):
                if i == 0:
                    self.Yout[i, :, 0] = y0_val_lst
                else:
                    self.Yout[i, :, 0] = self._r.integrate(self.tout[i])

                if self.nderiv > 0:
                    self.Yout[i, :, 1] = self._fo_odesys.dydt(
                        self.tout[i], self.Yout[i, :, 0], param_val_lst)
                if self.nderiv > 1:
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
                if self.nderiv > 0:
                    dyout.append(self._fo_odesys.dydt(
                        self._r.t, self._r.y, param_val_lst))
                if self.nderiv > 1:
                    ddyout.append(self._fo_odesys.d2ydt2(
                        self._r.t, self._r.y, param_val_lst))
                self._r.integrate(tend, step=True)
                yout.append(self._r.y)
                tout.append(self._r.t)
            warnings.resetwarnings()
            tout = np.array(tout)
            yout, dyout, ddyout = map(np.array, (yout, dyout, ddyout))
            if self.nderiv == 0:
                self.Yout = yout.reshape((yout.shape[0], yout.shape[1], 1))
            elif self.nderiv == 1:
                self.Yout = np.concatenate((yout[...,np.newaxis], dyout[...,np.newaxis]), axis=2)
            else:
                self.Yout = np.concatenate((yout[...,np.newaxis], dyout[...,np.newaxis],
                                           ddyout[...,np.newaxis]), axis=2)
            self.tout = tout

class SympyEvalr(object):
    """
    Specialized class to mimic the interface of IVP_Integrator
    class. It is used when an ODE has an analytic solution.

    Instances evaluates sympy expressions dependent on one variable and
    also substitutes parameters in expression using provided instance
    of FirstOrderODESystem
    The instances also evaluates the derivates in the independent variable
    up to requested order. (to facilitate interpolation)
    """

    default_dtype = np.float64

    def __init__(self, nderiv=0, dtype=None):
        """

        Arguments:
        - `exprs`: List of expressions to be evaluated
        - `indep_var_symb`: Sympy symbol of the independent variable
        - `params_by_symb`: Dictionary mapping sympy symbols of parameters to values
        - `nderiv`: set higher than 0 (default) in nderiv to also evaluate derivatives.
                   (output is sotred in self.Yout)
        """
        self.nderiv = nderiv
        if dtype == None: dtype = self.default_dtype
        self._dtype = dtype


    def configure(self, fo_odesys, param_vals):
        self._exprs = fo_odesys.solved_exprs
        self._indep_var_symb = fo_odesys.indepv
        self._params_by_symb = param_vals


    def eval_for_indep_array(self, arr, extra_params_by_symb):
        """
        Evaluate all expressions for values of indepedndent variable
        in array `arr` using self._params_symb and `extra_params_by_symb`
        for static substitution in sympy expressions in the list self._exprs
        """
        _Yout = np.empty((len(arr), len(self._exprs), self.nderiv+1), dtype=self._dtype)
        for expr_idx, expr in enumerate(self._exprs):
            for nderiv in range(self.nderiv+1):
                subs = self._params_by_symb
                subs.update(extra_params_by_symb)
                diff_expr = expr.diff(self._indep_var_symb, nderiv)
                print diff_expr
                for i, t in enumerate(arr):
                    subs.update({self._indep_var_symb: t})
                    _Yout[i, expr_idx, nderiv] = diff_expr.subs(subs)
        self.Yout = _Yout
