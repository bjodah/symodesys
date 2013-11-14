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
    h_init = None
    h_max  = 0.0 # inifinte

    _fo_odesys = None

    # Default highest order of derivatives to save for output
    nderiv = 1

    info = None

    def __init__(self, **kwargs):
        self.info = self.info or {}
        for key, val in kwargs.items():
            if hasattr(self, key):
                # convenient...
                setattr(self, key, val)

    def set_fo_odesys(self, fo_odesys):
        self._fo_odesys = fo_odesys


    def run(self, depv_init, params, indepv_init,
            indepv_end, N, check_jac_cond=False, **kwargs):
        """
        Run the numeric integration.

        Arguments:
        -`depv_init`: dictionary mapping dependent variables (string
          of its name or sympy.Symbol instance) to numeric values.
        -`params`: dictionary mapping parameters (string
          of its name or sympy.Symbol instance) to numeric values. Note
          that the parameters may be different from those of the
          original ODESys instance due to manipulation of the system
          inbetween (analytic treatment in particular)
        -`indepv_init`: numeric value of the initial value of the
          independent variable
        -`indepv_end`: the final value up to
          which integration in the independent variable is to be
          conducted.

        self._run assigns to self.tout and self.yout

        """
        # The C/Fortran/etc.. program knows nothing about dicts, provide values
        # as an array
        # Set (program) init vals to the (subset of non_anlytic) depv
        depv_init_arr = np.array(
            [depv_init[k] for k in self._fo_odesys.na_depv],
            dtype = np.float64)

        params_arr = np.array(
            [params[k] for k in self._fo_odesys.param_and_sol_symbs],
            dtype = np.float64)

        # Below is a somewhat ad-hoc sanity check of condition
        # of jacobian in the starting point, even though some
        # solvers might not need the jacobian at the starting point.
        # To work around this one would ideally use a variable transformation
        # and/or solving/estimating parts of the problem analytically
        if check_jac_cond:
            jac_cond = np.linalg.cond(self._fo_odesys.evaluate_na_jac(
                indepv_init, depv_init_arr, params_arr))#
            self.info['init_jac_cond'] = jac_cond
            if jac_cond*np.finfo(np.float64).eps > max(self.abstol, self.reltol):
                raise RuntimeError(("Unlikely that Jacboian with condition: {} "+\
                                   "will work with requested tolerances.").format(
                                       jcond))

        if self.h_init == None:
            self.h_init = 1e-9 # TODO: along the lines of:
            #   h_init=calc_h_init(depv_init, dydt, jac, abstol, reltol)
            # Also print a warning if h_init is below machine_epsilon**0.5
            # (or something) suggesting variable transformation.

        self._run(depv_init_arr, indepv_init, indepv_end, params_arr, N,
                  **kwargs)

    def clean(self):
        """
        Clean up temporary files.
        """
        pass


    def init_Yout_tout_for_fixed_step_size(self, indepv_init, indepv_end, N):
        dt = (indepv_end-indepv_init) / (N-1)
        NY = len(self._fo_odesys.na_depv)
        self.tout = np.asarray(np.linspace(indepv_init, indepv_end, N), dtype = self._dtype)
        # Handle other dtype for tout here? linspace doesn't support dtype arg..
        self.Yout = np.zeros((N, NY, self.nderiv+1), dtype = self._dtype)


class Mpmath_IVP_Integrator(IVP_Integrator):
    """
    Only for demonstration purposes - using mpmath has a severe performance
    penalty
    """

    def _run(self, depv_init_arr, indepv_init, indepv_end, params_arr, N):
        cb = lambda x, y: self._fo_odesys.evaluate_na_f(x, y, params_arr)
        self._num_y = sympy.mpmath.odefun(cb, indepv_init, depv_init_arr,
                                          tol = self.abstol)
        if N > 0:
            # Fixed stepsize
            self.init_Yout_tout_for_fixed_step_size(indepv_init, indepv_end, N)
            for i, t in enumerate(self.tout):
                if i == 0:
                    self.Yout[i, :, 0] = depv_init_arr
                else:
                    self.Yout[i, :, 0] = self._num_y(self.tout[i])

                if self.nderiv > 0:
                    self.Yout[i, :, 1] = self._fo_odesys.evaluate_na_f(
                        self.tout[i], self.Yout[i, :, 0], params_arr)
                if self.nderiv > 1:
                    self.Yout[i, :, 2] = self._fo_odesys.evaluate_na_d2ydt2(
                        self.tout[i], self.Yout[i, :, 0], params_arr)

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

    # def __init__(self, **kwargs):
    #     super(SciPy_IVP_Integrator, self).__init__(**kwargs)

    def set_fo_odesys(self, fo_odesys):
        super(SciPy_IVP_Integrator, self).set_fo_odesys(fo_odesys)
        from scipy.integrate import ode
        self._r = ode(self._fo_odesys.evaluate_na_f, self._fo_odesys.evaluate_na_jac)

    def _run(self, depv_init_arr, indepv_init, indepv_end, params_arr, N):
        self._r.set_initial_value(depv_init_arr, indepv_init)
        self._r.set_f_params(params_arr)
        self._r.set_jac_params(params_arr)
        if N > 0:
            # Fixed stepsize
            self._r.set_integrator('vode', method = 'bdf', with_jacobian = True)
            self.init_Yout_tout_for_fixed_step_size(indepv_init, indepv_end, N)
            for i, t in enumerate(self.tout):
                if i == 0:
                    self.Yout[i, :, 0] = depv_init_arr
                else:
                    self.Yout[i, :, 0] = self._r.integrate(self.tout[i])

                if self.nderiv > 0:
                    self.Yout[i, :, 1] = self._fo_odesys.evaluate_na_f(
                        self.tout[i], self.Yout[i, :, 0], params_arr)
                if self.nderiv > 1:
                    self.Yout[i, :, 2] = self._fo_odesys.evaluate_d2ydt2(
                        self.tout[i], self.Yout[i, :, 0], params_arr)
                assert self._r.successful()
        else:
            # Adaptive step size reporting
            # http://stackoverflow.com/questions/12926393/\
            #   using-adaptive-step-sizes-with-scipy-integrate-ode
            self._r.set_integrator('vode', method = 'bdf', with_jacobian = True, nsteps = 1)
            self._r._integrator.iwork[2] =- 1
            tout, yout, dyout, ddyout = [], [], [], []
            warnings.filterwarnings("ignore", category=UserWarning)
            keep_going = True
            while keep_going:
                keep_going = self._r.t < indepv_end
                if self.nderiv > 0:
                    dyout.append(self._fo_odesys.evaluate_na_f(
                        self._r.t, self._r.y, params_arr))
                if self.nderiv > 1:
                    ddyout.append(self._fo_odesys.evaluate_d2ydt2(
                        self._r.t, self._r.y, params_arr))
                self._r.integrate(indepv_end, step=True)
                yout.append(self._r.y)
                tout.append(self._r.t)
            warnings.resetwarnings()
            tout, yout, dyout, ddyout = map(np.array, (tout, yout, dyout, ddyout))
            if self.nderiv == 0:
                self.Yout = yout.reshape((yout.shape[0], yout.shape[1], 1))
            elif self.nderiv == 1:
                self.Yout = np.concatenate((yout[...,np.newaxis], dyout[...,np.newaxis]), axis=2)
            else:
                self.Yout = np.concatenate((yout[...,np.newaxis], dyout[...,np.newaxis],
                                           ddyout[...,np.newaxis]), axis=2)
            self.tout = tout
        self.info.update(self._fo_odesys.info)


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

    _dtype = np.float64

    nderiv = None

    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            if hasattr(self, key):
                # convenient...
                setattr(self, key, val)


    def set_fo_odesys(self, fo_odesys):
        self._fo_odesys = fo_odesys


    def eval_for_indep_array(self, indepv_arr, params):
        """
        Evaluate all expressions for values of indepedndent variable
        in array `indepv_arr` and using params (dict mapping Symbols to value)
        for static substitution in sympy expressions in the list
        self._fo_odesys.analytic_realtions
        """
        rels = self._fo_odesys.analytic_relations
        indepv = self._fo_odesys.indepv

        _Yout = np.empty((len(indepv_arr), len(rels), self.nderiv+1), dtype=self._dtype)
        for rel_idx, rel in enumerate(rels):
            for ideriv in range(self.nderiv+1):
                diff_expr = rel.diff(indepv, ideriv)
                for i, t in enumerate(indepv_arr):
                    params.update({indepv: t})
                    _Yout[i, rel_idx, ideriv] = diff_expr.subs(params)
        self.Yout = _Yout
