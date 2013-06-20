#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections import OrderedDict

import numpy as np
try:
    from cInterpol import PiecewisePolynomial
except ImportError:
    from scipy.interpolate import PiecewisePolynomial
import matplotlib.pyplot as plt
import sympy

from symodesys.helpers import cache, array_subs
from symodesys.integrator import Mpmath_IVP_Integrator, SympyEvalr
from symodesys.odesys import FirstOrderODESystem
from symodesys.transform import Transformer

# LONG-TERM FUTURE TODO: Support uncertainties as parameter inputs

def determine_const_val_for_init_val(expr, y0, indep_val,
                                     const_symb = sympy.symbols('C1')):
    """
    Helper function for IVP.recursive_analytic_reduction
    """
    y0_expr = expr.subs({indep_val: 0})
    const_val = sympy.solve(sympy.Eq(y0_expr, y0), const_symb)[0]
    return const_val


class IVP(object):
    """
    Initial Value Problem class

    The class abstracts away the change of initial values
    into parameters when a system of ODEs are reduced analytically
    I.e. the user may still update initial_values, even though in
    ``reality'' it is a parameter which is updated

    It also provides a front-end for the variable transformation
    routines and can automatically transforms numerical initial values
    and numerical results.

    Scaling of variables is provided in a similar manner.
    """

    # default_N is used if we integrate(..., N = 0, ...)
    # and all analytic sol. (no stepper is run)
    default_N = 100

    # Default highest order of derivatives to save for output
    nderiv = 1

    _dtype = np.float64

    def __init__(self, fo_odesys, init_vals, param_vals, t0,
                 integrator=None, analytic_evalr=None,
                 indepv_inv_trnsfm=None, depv_inv_trnsfm=None):
        """

        Arguments:
        - `fo_odesys`: First order ODE System
        - `init_vals`: Dictionary mapping dep. var symbols to vals at t0
        - `integrator`: IVP_Integrator class instance
        - `AnalyticEvalr`: Callback evaluating the analytically solved eq.
                            Defaults to SympyEvalr
        - `inv_indepv_trnfsm`: If provided output of numerical routines will be
           converted according to provided transform (usually provided by
           the method `use_internal_indepv_trnsfm`)
        - `inv_indepv_trnfsm`: If provided output of numerical routines will be
           converted according to provided transform (usually provided by
           the method `use_internal_depv_trnsfm`)

        """
        self._fo_odesys = fo_odesys
        self._old_fo_odesys = [] # Save old sys when solving
                                 # analytically
        self.init_vals = init_vals
        self.param_vals = param_vals
        self._indepv_init_val = t0
        self.integrator = integrator or Mpmath_IVP_Integrator()
        self.integrator.set_fo_odesys(self._fo_odesys)
        self.analytic_evalr = analytic_evalr or SympyEvalr()
        self.analytic_evalr.configure(self._fo_odesys, self.param_vals)
        self._indepv_inv_trnsfm = indepv_inv_trnsfm
        self._depv_inv_trnsfm = depv_inv_trnsfm

        # Init attributes for possible analytically solvable y's

        # TODO move _solved to ODESystem and handle accordingly..
        self._solved_init_val_symbs = {}


    def is_stiff(self, t0, tend, criteria=1e2):
        """
        Queries system using inintal values.
        """
        # Not working yet... (This is pseudo code mixed with to be fixed python)
        y0_val_lst = [self.init_vals[k] for k in self._fo_odesys.non_analytic_depv]
        param_val_lst = self._fo_odesys.param_val_lst(param_vals)
        return self._fo_odesys.is_stiff(t0, y0_val_lss, param_val_lst, criteria)


    @property
    def init_vals(self):
        return self._init_vals

    @init_vals.setter
    def init_vals(self, value):
        self._init_vals = self._fo_odesys.ensure_dictkeys_as_symbs(value)

    @property
    def param_vals(self):
        return self._param_vals

    @param_vals.setter
    def param_vals(self, value):
        self._param_vals = self._fo_odesys.ensure_dictkeys_as_symbs(value)

    def use_internal_depv_trnsfm(self, trnsfm, inv_trnsfm, strict=True):
        """
        Solve the system numerically for the transformed variables
        according to the provided arguments:
        -`trnsfm`: dict mapping new_depv to expr_in_old
        -`inv_trnsfm`: dict mapping old_depv to expression in new_depv
        -`strict`: make assertions that inv_trnsfm really is correct

        The user input and output of initial and resulting values
        for the dependent variables will be independent of the
        transformation (the values will be converted internally)
        """
        assert len(trnsfm) == len(self._fo_odesys.all_depv)
        assert len(trnsfm) == len(inv_trnsfm)

        if strict:
            for new_depv, expr_in_old in trnsfm.items():
                assert expr_in_old.subs(inv_trnsfm).simplify() == new_depv

        self._depv_trnsfm = trnsfm
        self._depv_inv_trnsfm = inv_trnsfm
        new_fo_odesys = self._fo_odesys.transform_depv(trnsfm, inv_trnsfm)
        new_init_vals = {}
        for new_depv, expr_in_old in trnsfm.items():
            new_init_vals[new_depv] = expr_in_old.subs(self.init_vals)
        return self.__class__(
            new_fo_odesys, new_init_vals, self.param_vals,
            self._indepv_init_val, self.integrator,
            self.analytic_evalr, self._indepv_inv_trnsfm, inv_trnsfm)


    def use_internal_indepv_trnsfm(self, trnsfm, inv_trnsfm, strict=True):
        """
        Returns a new class instance with this instance as model,
        but with independent variables transformed according
        to provided arguments
        -`trnsfm`: A tuple of (new_indep_symb, expr_in_old_indep)
        -`inv_trnsfm`: A tuple of (old_indep_symb, expr_in_new_indep)
        -`strict`: make assertions that inv_trnsfm really is correct
        """

        # look into https://groups.google.com/forum/?fromgroups#!topic/sympy/DNm2SpOdNd0
        if False: #strict:
            new_indep_symb, expr_in_old_indep = trnsfm
            assert expr_in_old_indep.subs(inv_trnsfm).simplify() == new_indep_symb

        self._indepv_trnsfm = trnsfm
        self._indepv_inv_trnsfm = inv_trnsfm
        new_fo_odesys = self._fo_odesys.transform_indepv(
            trnsfm, inv_trsfm)
        new_indepv_init_val = self._indepv_init_val
        return self.__class__(
            new_fo_odesys, self.init_vals, self.param_vals, new_indepv_init_val,
            self.integrator, self.analytic_evalr, inv_trnsfm, self._indepv_inv_trnsfm)

    def mk_depv_symbol(self, s):
        """
        Convenience function for generating sympy.Function symbols for
        use in `use_internal_depv_trnsfm`
        """
        return self._fo_odesys.mk_func(s)


    def mk_init_val_symb(self, y):
        new_symb = sympy.symbols(y.func.__name__ + '_init')
        assert not new_symb in self._fo_odesys.known_symbs
        return new_symb


    def recursive_analytic_reduction(self):
        """
        Attempts to solve some y's analytically

        TODO: recreate possible 2nd order ODE and check solvability
        """
        new_init_val_symbs = []
        self._fo_odesys.recursive_analytic_auto_sol()
        for yi, (expr, sol_symbs) in self._fo_odesys._solved.iteritems():
            if yi in self._solved_init_val_symbs: continue
            assert len(sol_symbs) == 1 # Only one new constant per equation
            sol_symb = sol_symbs.copy().pop()
            init_val_symb = self.mk_init_val_symb(yi)
            sol_init_val = determine_const_val_for_init_val(
                expr, init_val_symb, self._fo_odesys.indepv, sol_symb)
            self._fo_odesys.subs({sol_symb: sol_init_val})
            self._solved_init_val_symbs[yi] = init_val_symb
            new_init_val_symbs.append(init_val_symb)
        if len(new_init_val_symbs) > 0:
            self.analytic_evalr.configure(self._fo_odesys,
                                          self.param_vals)
        return new_init_val_symbs


    def integrate(self, tend, N=0, h=None, nderiv=None, **kwargs):
        """
        Integrates the non-analytic odesystem and evaluates the
        analytic functions for the dependent variables (if there
        are any).
        """
        assert float(tend) == tend
        self.indep_out.cache_clear()
        self._Yres.cache_clear()
        self.interpolators.cache_clear()
        self.trajectories.cache_clear()
        self.nderiv = nderiv or self.nderiv


        if len(self._solved_init_val_symbs) < len(self.init_vals):
            # If there are any non-analytic equations left
            y0 = {yi: self.init_vals[yi] for yi \
                  in self._fo_odesys.non_analytic_depv}
            self.integrator.run(
                y0, t0=self._indepv_init_val, tend=tend,
                param_vals=self.param_vals, N=N, h=h,
                nderiv=self.nderiv, **kwargs)
            #self.tout = self.integrator.tout
        else:
            # If all equations were solved analytically
            if N == 0: N = self.default_N
            self.integrator.tout = np.linspace(self._indepv_init_val, tend, N)

        if len(self._solved_init_val_symbs) > 0:
            self.analytic_evalr.nderiv = self.nderiv
            self.analytic_evalr.eval_for_indep_array(
                self.indep_out(), {
                    self._solved_init_val_symbs[yi]: self.init_vals[yi]\
                    for yi in self._fo_odesys.analytic_depv}
                )

    @cache
    def indep_out(self):
        """
        Handles variable transformation of numerical
        data corresponding to the independent variable
        """
        tout = self.integrator.tout
        if self._indepv_inv_trnsfm:
            new_tout = tout.copy()
            raise NotImplementedError
        else:
            return tout

    @cache
    def trajectories(self):
        """
        Returns an OrderedDict instance of:
        depv: nt×nderiv
        """
        Yres = self._Yres()
        tout = self.integrator.tout
        nt, ndepv, ndatapp = Yres.shape # time, dependent variables, data per point (nderiv+1)
        indepv = self._fo_odesys.indepv
        if self._depv_inv_trnsfm:

            # TODO: Handle self._indepv_inv_trnsfm

            # Ok, wee need to transform Yres from numerical integration
            # back to original variables
            od = OrderedDict()
            deriv_data = {depv.diff(indepv, j): Yres[:,i,j] for j in range(ndatapp) \
                           for i, depv in enumerate(self._fo_odesys.all_depv)}
            deriv_data[indepv] = tout

            exprs = []
            ori_derivs = []
            for ori_depv, expr_in_cur in self._depv_inv_trnsfm.items():
                for j in range(ndatapp):
                    der_expr = expr_in_cur.diff(indepv, j)
                    exprs.append(der_expr)
                    ori_derivs.append(ori_depv.diff(indepv, j))

            inp, yres_data = deriv_data.keys(), deriv_data.values()
            tmfr = Transformer(exprs, inp)
            tmfr_data = tmfr(*yres_data)
            for ori_depv in self._depv_inv_trnsfm.keys():
                idxs = []
                for i in range(ndatapp):
                    idxs.append(ori_derivs.index(ori_depv.diff(indepv, i)))
                od[ori_depv] = tmfr_data[:,idxs]
            return od
        else:
            return OrderedDict(zip(self._fo_odesys.all_depv,
                                   [Yres[:,i,:] for i in range(ndepv)]))


    @cache
    def _Yres(self):
        """
        Unified the output of the numerical and analyitc results.
        first axis: independent variable value
        second axis: dependent variable index (_fo_odesys.all_depv)
        third axis: 0-th, 1st, ... derivatives
        """
        Yres = np.empty((len(self.indep_out()), len(self._fo_odesys.all_depv),
                          self.nderiv+1), self._dtype)
        for i, yi in enumerate(self._fo_odesys.all_depv):
            if yi in self._fo_odesys.analytic_depv:
                Yres[:, i, :] = self.analytic_evalr.Yout[
                    :, self._fo_odesys.analytic_depv.index(yi),:]
            else:
                if len(self._fo_odesys.non_analytic_depv) > 0:
                    Yres[:, i, :] = self.integrator.Yout[
                        :, self._fo_odesys.non_analytic_depv.index(yi),:]
        return Yres

    @property
    def all_depv(self):
        """
        Resolves current depv (in case of internal variables transformation in use)
        """
        if self._depv_inv_trnsfm:
            return self._depv_inv_trnsfm.keys()
        else:
            return self._fo_odesys.all_depv

    @cache
    def interpolators(self):
        return OrderedDict([(k, PiecewisePolynomial(
            self.indep_out(), self.trajectories()[k])) for k,v \
                            in self.trajectories().items()])


    def get_interpolated(self, t, depvs=None):
        if depvs == None: depvs = self.all_depv
        return np.array([self.interpolators()[depv](t) for depv in depvs])


    def get_depv_from_token(self, depvn):
        """
        Like __getitem__ of FirstOrderODESystem, but intercepts
        variable use of internal variable transformation.
        """
        if self._depv_inv_trnsfm:
            candidate = self.mk_depv_symbol(depvn)
            if candidate in self._depv_inv_trnsfm.keys():
                return candidate
            else:
                raise KeyError('{} (created from {}) not found in original depv'.format(
                    candidate, depvn))
        else:
            return self._fo_odesys[depvn]


    def get_index_of_depv(self, depvn):
        return self.all_depv.index(self.get_depv_from_token(depvn))

    def plot(self, depvs=None, interpolate=True, datapoints=False,
             show=False, skip_helpers=True, usetex=False, texnames=None,
             ax=None):
        """
        Rudimentary plotting utility for quick inspection of solutions
        TODO: move to symodesys.convenience ?
        Arguments:
        - `depvs`: A sequence of depv to be plotted, plots all if it is None (default)
        """
        import matplotlib.pyplot as plt
        if usetex:
            from matplotlib import rc
            rc('text', usetex=True)

        if depvs == None:
            depvs = self.all_depv
            if skip_helpers:
                # Don't plot helper functions used in reduction of order of ode system
                for hlpr in self._fo_odesys.frst_red_hlprs:
                    depvs.pop(depvs.index(hlpr[2]))
        if interpolate:
            ipx = np.linspace(self.indep_out()[0], self.indep_out()[-1], 1000)
            ipy = self.get_interpolated(ipx, depvs)
        ls = ['-', '--', ':']
        c = 'k b r g m'.split()
        m = 'o s ^ * d p h'.split()
        ax = ax or plt.subplot(111)

        for i, depv in enumerate(depvs):
            mi  = m[i % len(m)]
            lsi = ls[i % len(ls)]
            ci  = c[i % len(c)]
            # Assign label
            lbl = str(depv)
            if usetex and texnames != None:
                lbl = texnames[lbl]
            if interpolate:
                ax.plot(ipx, ipy[i,:], label = lbl,
                         marker = 'None', ls = lsi, color = ci)
                lsi = 'None'
            if datapoints:
                ax.plot(self.indep_out(), self.trajectories()[depv][:, 0], label = lbl,
                         marker = mi, ls = lsi, color = ci)


        # Put a legend to the right of the current axis
        if show:
            # Shrink box by 20%
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            plt.show()
        return ax

    def __enter__(self): return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """
        Clean
        """
        # implement a clean method in the subclass if
        # you wish to use the context manager.
        self.clean()

    def clean(self):
        self.integrator.clean()
