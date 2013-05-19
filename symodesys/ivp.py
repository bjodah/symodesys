from collections import OrderedDict

import numpy as np
try:
    from cinterpol import PiecewisePolynomial
except ImportError:
    from scipy.interpolate import PiecewisePolynomial
import matplotlib.pyplot as plt
import sympy

from symodesys.helpers import SympyEvalr, cache
from symodesys.integrator import Mpmath_IVP_Integrator
from symodesys.odesys import FirstOrderODESystem

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

    _dtype = np.float64

    def __init__(self, fo_odesys, init_vals, param_vals, t0,
                 integrator = None, analytic_evalr = None):
        """

        Arguments:
        - `fo_odesys`: First order ODE System
        - `init_vals`: Dictionary mapping dep. var symbols to vals at t0
        - `integrator`: IVP_Integrator class instance
        - `AnalyticEvalr`: Callback evaluating the analytically solved eq.
                            Defaults to SympyEvalr
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

        # Init attributes for possible analytically solvable y's

        # TODO move _solved to ODESystem and handle accordingly..
        self._solved_init_val_symbs = {}
        self._depv_trnsfm = None
        self._depv_inv_trnsfm = None
        self._indepv_trnsfm = None
        self._indepv_inv_trnsfm = None


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


    def use_internal_depv_trnsfm(self, trnsfm, inv_trnsfm):
        """
        Solve the system numerically for the transformed variables
        according to the provided arguments:
        -`trnsfm`: dict mapping new_depv to expr_in_old
        -`inv_trnsfm`: dict mapping old_depv to expression in new_depv

        The user input and output of initial and resulting values
        for the dependent variables will be agnostic of the
        transformation (the values will be converted internally)
        """
        self._depv_trnsfm = trnsfm
        self._depv_inv_trnsfm = inv_trnsfm
        new_fo_odesys = self._fo_odesys.transform_depv(trnsfm, inv_trnsfm)
        new_init_vals = {}
        for new_depv, expr_in_old in trnsfm.items():
            new_init_vals[new_depv] = expr_in_old.subs(self.init_vals)
        new_instance =  self.__class__(new_fo_odesys, new_init_vals,
                              self.param_vals, self._indepv_init_val,
                              self.integrator, self.analytic_evalr)
        new_instance._inv_depv_trnsf = inv_trnsfm
        return new_instance


    def use_internal_indepv_trnsfm(self, trnsfm, inv_trnsfm):
        """
        Returns a new class instance with this instance as model,
        but with independent variables transformed according
        to provided arguments
        -`trnsfm`: A tuple of (new_indep_symb, expr_in_old_indep)
        -`inv_trnsfm`: A tuple of (old_indep_symb, expr_in_new_indep)
        """
        self._indepv_trnsfm = trnsfm
        self._indepv_inv_trnsfm = inv_trnsfm
        new_fo_odesys = self._fo_odesys.transform_indepv(
            trnsfm, inv_trsfm)
        new_indepv_init_val = self._indepv_init_val
        new_instance = self.__class__(new_fo_odesys, self.init_vals,
                              self.param_vals, new_indepv_init_val,
                              self.integrator, self.analytic_evalr)
        new_instance._inv_indepv_trnsfm = inv_trnsfm
        return new_instance

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
        if len(new_init_val_symbs) > 0:
            self.analytic_evalr.configure(self._fo_odesys,
                                          self.param_vals)


    def integrate(self, tend, N = 0, h = None, order = 1):
        """
        Integrates the non-analytic odesystem and evaluates the
        analytic functions for the dependent variables (if there
        are any).
        """
        assert float(tend) == tend
        self.interpolators.cache_clear()
        self.trajectories.cache_clear()

        if len(self._solved_init_val_symbs) < len(self.init_vals):
            # If there are any non-analytic equations left
            y0 = {yi: self.init_vals[yi] for yi \
                  in self._fo_odesys.non_analytic_depv}
            self.integrator.run(y0,
                t0 = self._indepv_init_val, tend = tend,
                param_vals = self.param_vals,
                N = N, h = h, order = order)
            #self.tout = self.integrator.tout
        else:
            if N == 0: N = self.default_N
            #self.tout = np.linspace(self._indepv_init_val, tend, N)

        if len(self._solved_init_val_symbs) > 0:
            self._analytic_evalr.eval_for_indep_array(
                self.tout, {
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
        if hasattr(self, '_inv_indepv_trnsfm'):
            new_tout = tout.copy()
        else:
            return tout


    @cache
    def trajectories(self):
        """
        Handles variable transformation of numerical
        data corresponding to the dependent variables
        """
        Yres = self.Yres
        if hasattr(self, '_inv_depv_trnsfm'):
            # do the inverse transform of the dependent variables.
            new_Yres=Yres.copy()
            Yres_dict={}
            for i, cur_depv in enumerate(self._fo_odesys.all_depv):
                for j in rang(Yres.shape[2]):
                    # j loops over ith deriv
                    Yres_dict[cur_depv.diff(self._fo_odesys.indepv, i)] \
                        = Yres[:, i, j]

            for i, (ori_depv, expr_in_cur) in enumerate(
                    self._inv_depv_trnsfm.items()):
                derivs = []
                for j in range(Yres.shape[2]):
                    # j loops over ith deriv
                    der_expr = expr_in_cur.diff(
                        self._fo_odesys.indepv)
                    for k in range(Yres.shape[0]):
                        # ouch, this will be slow
                        new_Yres[k,i,j] = der_expr.subs(Yres_dict)

            return new_Yres
        else:
            return Yres


    @cache
    def _Yres(self):
        """
        Unified the output of the numerical and analyitc results.
        first axis: independent variable value
        second axis: dependent variable index (_fo_odesys.all_depv)
        third axis: 0-th, 1st, ... derivatives
        """
        if not hasattr(self, 'tout'): return None
        Yres = np.empty((len(self.tout), len(self._fo_odesys.all_depv),
                          self.integrator.Yout.shape[2]), self._dtype)
        for i, yi in enumerate(self._fo_odesys.all_depv):
            if yi in self._fo_odesys.analytic_depv:
                Yres[:, i, :] = self._analytic_evalr.Yout[
                    :, self._fo_odesys.analytic_depv.index(yi),:]
            else:
                Yres[:, i, :] = self.integrator.Yout[
                    :, self._fo_odesys.non_analytic_depv.index(yi),:]
        return Yres

    def depv_indices(self):
        return range(self.trajectories().shape[1])


    @cache
    def interpolators(self):
        intrpltrs = []
        for i in self.depv_indices():
            intrpltrs.append(PiecewisePolynomial(self.tout, self.trajectories()[:,i,:]))
        return intrpltrs

    def get_interpolated(self, t):
        return np.array([self.interpolators()[i](t) for i in self.depv_indices()])

    def get_index_of_depv(self, depvn):
        return self._fo_odesys.all_depv.index(self._fo_odesys[depvn])

    def plot(self, indices = None, interpolate = True, datapoints=False,
             show = False, skip_helpers = True, usetex=False, texnames=None):
        """
        Rudimentary plotting utility for quick inspection of solutions
        TODO: move this from here,  make more general to accept mixed ODE sol +
        analytic y curves
        """
        import matplotlib.pyplot as plt
        if usetex:
            from matplotlib import rc
            rc('text', usetex=True)

        if indices == None:
            indices = self.depv_indices()
            if skip_helpers:
                # Don't plot helper functions used in reduction of order of ode system
                for hlpr in self._fo_odesys.frst_red_hlprs:
                    indices.pop(self._fo_odesys.all_depv.index(hlpr[2]))
        if interpolate:
            ipx = np.linspace(self.tout[0], self.tout[-1], 1000)
            ipy = self.get_interpolated(ipx)
        ls = ['-', '--', ':']
        c = 'k b r g m'.split()
        m = 'o s ^ * d p h'.split()
        for i in indices:
            mi  = m[i % len(m)]
            lsi = ls[i % len(ls)]
            ci  = c[i % len(c)]
            lbl = str(self._fo_odesys.all_depv[i])
            if usetex and texnames != None:
                lbl = texnames[lbl]
            print(lbl) ###
            if interpolate:
                plt.plot(ipx, ipy[i,:], label = lbl + ' (interpol.)',
                         marker = 'None', ls = lsi, color = ci)
                lsi = 'None'
            if datapoints:
                plt.plot(self.tout, self.trajectories()[:, i, 0], label = lbl,
                         marker = mi, ls = lsi, color = ci)
            plt.plot()
        plt.legend()
        if show: plt.show()

def plot_numeric_vs_analytic(ODESys, y0, params, tend, t0=0.0, N=0):
    """
    Convenience function for instantiating ODESys class and assigning
    it to an associated IVP instance, run the integration and in the case
    of ODESys having 'analytic_sol', plot the absolute and relative errors
    made during the integration.
    """
    odesys = ODESys()
    ivp = IVP(odesys, y0, params, t0)
    ivp.integrate(tend, N)
    t, y = ivp.tout, ivp.trajectories()[:,:,0]

    plt.subplot(311)
    ivp.plot(interpolate = True, show = False)

    for i, (k, cb) in enumerate(odesys.analytic_sol.items()):
        analytic = cb(odesys, t, y0, params)
        plt.subplot(312)
        plt.plot(t, (y[:, i] - analytic) / ivp.integrator.abstol,
                 label = k+': abserr / abstol')
        plt.legend()
        plt.subplot(313)
        plt.plot(t, (y[:, i] - analytic) / analytic / ivp.integrator.reltol,
                 label = k+': relerr / reltol')
        plt.legend()

    plt.show()
