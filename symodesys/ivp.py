from collections import OrderedDict

import numpy as np
try:
    import cinterpol
except ImportError:
    cinterpol = False
from scipy.interpolate import PiecewisePolynomial
import matplotlib.pyplot as plt
import sympy

from symodesys.helpers import SympyEvalr, cache
from symodesys.integrator import SciPy_IVP_Integrator
from symodesys.firstorder import FirstOrderODESystem

# FUTURE: Support uncertainties as parameter inputs

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
    """

    # TODO: instead of self.yout, self.dyout, self.ddyout, ...
    #       implement self.Yout (3 dimensional, where
    #          first dimension is order)

    # used if integrate(..., N = 0, ...) and all analytic sol.
    default_N = 100

    _dtype = np.float64

    def __init__(self, fo_odesys, init_vals, param_vals_by_symb, t0,
                 Integrator = SciPy_IVP_Integrator,
                 AnalyticEvalr = SympyEvalr):
        """

        Arguments:
        - `fo_odesys`: First order ODE System
        - `init_vals`: Dictionary mapping dep. var symbols to vals at t0
        - `Integrator`: IVP_Integrator class
        - `AnalyticEvalr`: Callback evaluating the analytically solved eq.
                            Defaults to SympyEvalr
        """
        self._fo_odesys = fo_odesys
        self._old_fo_odesys = [] # Save old sys when solving analytically
        self._init_vals = init_vals
        self._param_vals_by_symb = param_vals_by_symb
        self._indepv_init_val = t0
        self._Integrator = Integrator
        self._AnalyticEvalr = AnalyticEvalr

        # Init attributes for possible analytically solvable y's

        # TODO move _solved to ODESystem and handle accordingly..
        self._solved_init_val_symbs = {}
        self._depv_trnsfm = None
        self._depv_inv_trnsfm = None
        self._indepv_trnsfm = None
        self._indepv_inv_trnsfm = None


    def use_internal_depv_trnsfm(self, trnsfm, inv_trnsfm):
        self._depv_trnsfm = trnsfm
        self._depv_inv_trnsfm = inv_trnsfm
        self._fo_odesys.transform_depv(trnsfm, inv_trnsfm)

    def use_internal_indepv_trnsfm(self, trnsfm, inv_trnsfm):
        self._indepv_trnsfm = trnsfm
        self._indepv_inv_trnsfm = inv_trnsfm
        self._fo_odesys.transform_indepv(trnsfm, inv_trsfm)

    def mk_init_val_symb(self, y):
        new_symb = sympy.symbols(y.func.__name__ + '_init')
        assert not new_symb in self._fo_odesys.known_symbs
        return new_symb


    def recursive_analytic_reduction(self):
        """
        Attempts to solve some y's analytically

        TODO: recreate possible 2nd order ODE and check solvability
        """
        new_init_val_symbs = [] # For info reporting
        self._fo_odesys.recursive_analytic_auto_sol()
        for yi, (expr, sol_symbs) in self._fo_odesys.solved.iteritems():
            if yi in self._solved_init_val_symbs: continue
            assert len(sol_symbs) == 1
            sol_symb = sol_symbs.copy().pop()
            init_val_symb = self.mk_init_val_symb(yi)
            sol_init_val = determine_const_val_for_init_val(
                expr, init_val_symb, self._fo_odesys.indepv, sol_symb)
            self._fo_odesys.subs({sol_symb: sol_init_val})
            self._solved_init_val_symbs[yi] = init_val_symb
            new_init_val_symbs.append(init_val_symb)
            return new_init_val_symbs


    def integrate(self, tend, N = 0, h = None, order = 2):
        """
        Integrates the non-analytic odesystem and evaluates the
        analytic functions for the dependent variables (if there are any).
        """
        if len(self._solved_init_val_symbs) > 0:
            # If we have solved parts analytically
            self._analytic_evalr = self._AnalyticEvalr(
                self._fo_odesys.solved_exprs,
                self._fo_odesys.indepv,
                self._param_vals_by_symb, order = order)

        print len(self._solved_init_val_symbs)
        print len(self._init_vals)
        if len(self._solved_init_val_symbs) < len(self._init_vals):
            # If there are any non-analytic equations left
            self._integrator = self._Integrator(self._fo_odesys)
            self._integrator.integrate(
                {yi: self._init_vals[yi] for yi \
                 in self._fo_odesys.non_analytic_depv},
                t0 = self._indepv_init_val, tend = tend,
                param_vals_by_symb = self._param_vals_by_symb,
                N = N, h = h, order = order)
            self.tout = self._integrator.tout
        else:
            if N == 0: N = self.default_N
            self.tout = np.linspace(self._indepv_init_val, tend, N)

        if len(self._solved_init_val_symbs) > 0:
            self._analytic_evalr.eval_for_indep_array(
                self.tout, {
                    self._solved_init_val_symbs[yi]: self._init_vals[yi]\
                    for yi in self._fo_odesys.analytic_depv}
                )

    @property
    def trajectory(self):
        """
        Handles inv_trnsfm
        """
        pass

    @property
    def yout(self):
        """
        Unified the output of the numerical and analyitc results.
        """
        _yout = np.empty((len(self.tout), len(self._fo_odesys.all_depv)),
                         self._dtype)
        for i, yi in enumerate(self._fo_odesys.all_depv):
            if yi in self._fo_odesys.analytic_depv:
                _yout[:, i] = self._analytic_evalr.yout[
                    :, self._fo_odesys.analytic_depv.index(yi)]
            else:
                _yout[:, i] = self._integrator.yout[
                    :, self._fo_odesys.non_analytic_depv.index(yi)]
        return _yout

    def Dy(self):
        _out = np.empty((len(self.tout), len(self._fo_odesys.all_depv)),
                        self._dtype)
        if len(self._fo_odesys.non_analytic_depv) > 0:
            num_dyout = np.array(
                [self._fo_odesys.dydt(
                    t, self._integrator.yout[i,:],
                    self._fo_odesys.param_val_lst(
                        self._param_vals_by_symb)
                    ) for\
                 (i,), t in np.ndenumerate(self.tout)])
        for i, yi in enumerate(self._fo_odesys.all_depv):
            if yi in self._fo_odesys.analytic_depv:
                _out[:, i] = self._analytic_evalr.dyout[
                    :, self._fo_odesys.analytic_depv.index(yi)]
            else:
                # if len(self._fo_odesys.all_depv) > 0:
                _out[:, i] = num_dyout[
                    :, self._fo_odesys.non_analytic_depv.index(yi)]
        return _out


    def DDy(self):
        nt = self.tout.shape[0]
        _out = np.empty((nt, len(self._fo_odesys.all_depv)), self._dtype)
        if len(self._fo_odesys.non_analytic_depv) > 0:
            num_ddyout = np.array(
                [self._fo_odesys.d2ydt2(
                    t, self.yout[i,:], self._fo_odesys.param_val_lst(
                        self._param_vals_by_symb)) for\
                 (i,), t in np.ndenumerate(self.tout)])
        for i, yi in enumerate(self._fo_odesys.all_depv):
            if yi in self._fo_odesys.analytic_depv:
                _out[:, i] = self._analytic_evalr.ddyout[
                    :, self._fo_odesys.analytic_depv.index(yi)]
            else:
                #if len(self._fo_odesys.all_depv) > 0:
                _out[:, i] = num_ddyout[
                    :, self._fo_odesys.non_analytic_depv.index(yi)]
        return _out

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
                     in self._fo_odesys.all_depv])[str(symb)][:, 0]

    def plot(self, indices = None, interpolate = False, show = False, skip_helpers = True):
        """
        Rudimentary plotting utility for quick inspection of solutions
        TODO: move this from here,  make more general to accept mixed ODE sol +
        analytic y curves
        """
        if indices == None:
            indices = range(self.yout.shape[1])
            if skip_helpers:
                # Don't plot helper functions used in reduction of order of ode system
                for hlpr in self._fo_odesys.frst_red_hlprs:
                    indices.pop(self._fo_odesys.all_depv.index(hlpr[2]))
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
            lbl = str(self._fo_odesys.all_depv[i])
            if interpolate:
                plt.plot(ipx, ipy[:, i], label = lbl + ' (interpol.)',
                         marker = 'None', ls = lsi, color = ci)
                lsi = 'None'
            plt.plot(self.tout, self.yout[:, i], label = lbl,
                     marker = mi, ls = lsi, color = ci)
            plt.plot()
        plt.legend()
        if show: plt.show()

