#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import, unicode_literals

from collections import OrderedDict
from functools import partial

import numpy as np
try:
    from cInterpol import PiecewisePolynomial
except ImportError:
    from scipy.interpolate import PiecewisePolynomial
import matplotlib.pyplot as plt
import sympy

from symodesys.helpers import cache
from symodesys.integrator import Mpmath_IVP_Integrator, SympyEvalr
from symodesys.odesys import FirstOrderODESystem
from symvarsub import NumTransformer


def determine_const_val_for_init_val(expr, depv, const_symb,
                                     indepv, indepv_init,
                                     init_symb_factory
                                     ):
    """
    Helper function for IVP.recursive_analytic_reduction
    """
    depv_init = init_symb_factory(depv)
    new_expr = expr.subs({indepv: indepv_init,
                          depv: depv_init})
    const_symb_in_init = sympy.solve(sympy.Eq(new_expr, depv_init), const_symb)[0]
    return const_symb_in_init, depv_init


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


    _dtype = np.float64

    _indepv_init_symb = None # used in analytic solution of IVP

    def __init__(self, fo_odesys, depv_init, params, indepv_init,
                 integrator=None, analytic_evalr=None,
                 indepv_inv_trnsfm=None, depv_inv_trnsfm=None,
                 logger=None):
        """

        Arguments:
        - `fo_odesys`: First order ODE System
        - `depv_init`: Dictionary mapping dep. var symbols to vals at t0
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
        self.fo_odesys = fo_odesys
        self._old_fo_odesys = [] # Save old sys when solving
                                 # analytically
        self.depv_init = depv_init
        self.params = params
        self._indepv_init_val = indepv_init
        self.integrator = integrator or Mpmath_IVP_Integrator()
        self.analytic_evalr = analytic_evalr or SympyEvalr(
            nderiv=self.integrator.nderiv)
        self._indepv_inv_trnsfm = indepv_inv_trnsfm
        self._depv_inv_trnsfm = depv_inv_trnsfm
        self.logger = logger


    def check_if_stiff(self, t0, tend, criteria=1e2):
        """
        Queries system using inintal values.
        """
        # Not working yet... (This is pseudo code mixed with to be fixed python)
        y0_val_lst = [self.depv_init[k] for k in self.fo_odesys.na_depv]
        param_val_lst = self.fo_odesys.param_val_lst(params)
        return self.fo_odesys.is_stiff(t0, y0_val_lss, param_val_lst, criteria)


    @property
    def depv_init(self):
        return self._depv_init

    @depv_init.setter
    def depv_init(self, value):
        self._depv_init = self.fo_odesys.ensure_dictkeys_as_symbs(value)

    @property
    def params(self):
        return self._params

    @params.setter
    def params(self, value):
        self._params = self.fo_odesys.ensure_dictkeys_as_symbs(value)

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
        assert len(trnsfm) == len(self.fo_odesys.all_depv)
        assert len(trnsfm) == len(inv_trnsfm)

        if strict:
            for new_depv, expr_in_old in trnsfm.items():
                assert expr_in_old.subs(inv_trnsfm).simplify() == new_depv

        self._depv_trnsfm = trnsfm
        self._depv_inv_trnsfm = inv_trnsfm
        new_fo_odesys = self.fo_odesys.transform_depv(trnsfm, inv_trnsfm)
        new_depv_init = {}
        for new_depv, expr_in_old in trnsfm.items():
            new_depv_init[new_depv] = expr_in_old.subs(self.depv_init)
        return self.__class__(
            new_fo_odesys, new_depv_init, self.params,
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
        new_fo_odesys = self.fo_odesys.transform_indepv(
            trnsfm, inv_trsfm)
        new_indepv_init_val = self._indepv_init_val
        return self.__class__(
            new_fo_odesys, self.depv_init, self.params, new_indepv_init_val,
            self.integrator, self.analytic_evalr, inv_trnsfm, self._indepv_inv_trnsfm)


    def mk_init_val_symb(self, y):
        new_symb = self.fo_odesys.mk_symb(y.func.__name__ + '_init')
        assert not new_symb in self.fo_odesys.known_symbs
        return new_symb


    @property
    def indepv_init_symb(self):
        if not self._indepv_init_symb:
            self._indepv_init_symb = self.fo_odesys.mk_symb(
                self.fo_odesys.indepv.name+'_init')
            self.fo_odesys.param_symbs.append(self._indepv_init_symb)
        return self._indepv_init_symb

    def recursive_analytic_reduction(self, complexity=0):
        """
        Attempts to solve some y's analytically

        TODO: recreate possible 2nd order ODE and check solvability
        """
        nsol = self.fo_odesys.recursive_analytic_auto_sol(
            complexity,
            cb=partial(determine_const_val_for_init_val,
                       indepv=self.fo_odesys.indepv,
                       indepv_init=self.indepv_init_symb,
                       init_symb_factory=self.mk_init_val_symb),
            logger=self.logger)

        if nsol > 0:
            self.analytic_evalr.set_fo_odesys(self.fo_odesys)
            if self.logger: self.logger.info('Done solving!')
        else:
            if self.logger: self.logger.info(
                    ('Unable to make any analytic reductions at'+\
                     ' complexity={}').format(complexity))


    def _clear_caches(self):
        self.indepv_out.cache_clear()
        self._Yres.cache_clear()
        self.interpolators.cache_clear()
        self.trajectories.cache_clear()


    def integrate(self, indepv_end, N, **kwargs):
        """
        Integrates the non-analytic odesystem and evaluates the
        analytic functions for the dependent variables (if there
        are any).
        """
        # sanity check
        assert float(indepv_end) == indepv_end

        self._clear_caches()

        # We might have introduced some new parameters,
        # make sure to include them in `params_complete`
        params_complete = self.params
        for depv, (expr, sol_symb) in self.fo_odesys._solved.items():
            params_complete.update({sol_symb[0]: self.depv_init[depv]})
        if self._indepv_init_symb:
            params_complete.update(
                {self._indepv_init_symb: self._indepv_init_val})

        if len(self.fo_odesys._solved) < len(self.depv_init):
            # If there are any non-analytic equations left
            self.integrator.set_fo_odesys(self.fo_odesys)
            self.integrator.run(
                self.depv_init, indepv_init=self._indepv_init_val,
                indepv_end=indepv_end, params=params_complete, N=N, **kwargs)
        else:
            # If all equations were solved analytically
            self.integrator.tout = np.linspace(self._indepv_init_val,
                                               indepv_end, N)

        if len(self.fo_odesys._solved) > 0:
            self.analytic_evalr.eval_for_indep_array(
                self.integrator.tout, params_complete
                )
        if self.logger: self.logger.info('Detailed info about integration:\n{}'.format(
                self.integrator.info))

    @cache
    def indepv_out(self):
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
        indepv = self.fo_odesys.indepv
        if self._depv_inv_trnsfm:

            # TODO: Handle self._indepv_inv_trnsfm

            # Ok, wee need to transform Yres from numerical integration
            # back to original variables
            od = OrderedDict()
            deriv_data = {depv.diff(indepv, j): Yres[:,i,j] for j in range(ndatapp) \
                           for i, depv in enumerate(self.fo_odesys.all_depv)}
            deriv_data[indepv] = tout

            exprs = []
            ori_derivs = []
            for ori_depv, expr_in_cur in self._depv_inv_trnsfm.items():
                for j in range(ndatapp):
                    der_expr = expr_in_cur.diff(indepv, j)
                    exprs.append(der_expr)
                    ori_derivs.append(ori_depv.diff(indepv, j))

            inp, yres_data = deriv_data.keys(), deriv_data.values()
            tmfr = NumTransformer(exprs, inp, save_temp=True) ###
            tmfr_data = tmfr(*yres_data)
            for ori_depv in self._depv_inv_trnsfm.keys():
                idxs = []
                for i in range(ndatapp):
                    idxs.append(ori_derivs.index(ori_depv.diff(indepv, i)))
                od[ori_depv] = tmfr_data[:,idxs]
            return od
        else:
            return OrderedDict(zip(self.fo_odesys.all_depv,
                                   [Yres[:,i,:] for i in range(ndepv)]))


    @cache
    def _Yres(self):
        """
        Unified the output of the numerical and analyitc results.
        first axis: independent variable value
        second axis: dependent variable index (fo_odesys.all_depv)
        third axis: 0-th, 1st, ... derivatives
        """
        Yres = np.empty((len(self.indepv_out()), len(self.fo_odesys.all_depv),
                          self.integrator.nderiv+1), self._dtype)
        for i, yi in enumerate(self.fo_odesys.all_depv):
            if yi in self.fo_odesys.analytic_depv:
                Yres[:, i, :] = self.analytic_evalr.Yout[
                    :, self.fo_odesys.analytic_depv.index(yi),:]
            else:
                if len(self.fo_odesys.na_depv) > 0:
                    Yres[:, i, :] = self.integrator.Yout[
                        :, self.fo_odesys.na_depv.index(yi),:]
        return Yres


    @property
    def all_depv(self):
        """
        Resolves current depv (in case of internal variables transformation in use)
        """
        if self._depv_inv_trnsfm:
            return self._depv_inv_trnsfm.keys()
        else:
            return self.fo_odesys.all_depv


    @cache
    def interpolators(self):
        return OrderedDict([(k, PiecewisePolynomial(
            self.indepv_out(), self.trajectories()[k])) for k,v \
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
            candidate = self.fo_odesys.mk_depv(depvn)
            if candidate in self._depv_inv_trnsfm.keys():
                return candidate
            else:
                raise KeyError('{} (created from {}) not found in original depv'.format(
                    candidate, depvn))
        else:
            return self.fo_odesys[depvn]


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
                for hlpr in self.fo_odesys.frst_red_hlprs:
                    depvs.pop(depvs.index(hlpr[2]))
        if interpolate:
            ipx = np.linspace(self.indepv_out()[0], self.indepv_out()[-1], 1000)
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
                ax.plot(self.indepv_out(), self.trajectories()[depv][:, 0], label = lbl,
                         marker = mi, ls = lsi, color = ci)

        if hasattr(self.fo_odesys, 'title'):
            if self.fo_odesys.title: ax.title(self.fo_odesys.title)

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
