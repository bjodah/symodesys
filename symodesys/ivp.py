from collections import OrderedDict

import numpy as np
from scipy.interpolate import PiecewisePolynomial
import matplotlib.pyplot as plt
import sympy

from symodesys.helpers import SympyEvalr, cache
from symodesys.integrator import SciPy_IVP_Integrator

def determine_const_val_for_init_val(f_general_sol, f0, dep_val,
                                     const_symb = sympy.symbols('C1')):
    """
    Helper function for IVP.recursive_analytic_reduction
    """
    f0_expr = f_general_sol.rhs.subs({dep_val: 0})
    const_val = sympy.solve(sympy.Eq(f0_expr, f0), const_symb)[0]
    return sympy.Eq(f_general_sol.lhs,
                    f_general_sol.rhs.subs({const_symb: const_val}))


class IVP(object):
    """
    Initial Value Problem class

    The class abstracts away the change of initial values
    into parameters when a system of ODEs are reduced analytically
    I.e. the user may still update initial_values, even though in
    ``reality'' it is a parameter which is updated
    """

    default_N = 100 # used if integrate(..., N = 0, ...) and all analytic sol.

    _dtype = np.float64

    def __init__(self, fo_odesys, init_vals,
                 Integrator = SciPy_IVP_Integrator,
                 AnalyticEvalr = SympyEvalr):
        """

        Arguments:
        - `fo_odesys`: First order ODE System
        - `init_vals`: Dictionary mapping dep. var names to initial values
        - `Integrator`: IVP_Integrator class
        - `AnalyticEvalr`: Callback evaluating the analytically solved eq.
                            Defaults to SympyEvalr
        """
        self._fo_odesys = fo_odesys
        self._init_vals = init_vals
        self._Integrator = Integrator
        self._AnalyticEvalr = AnalyticEvalr

        # Save the original dependent variables symbols:
        self._ori_dep_var_func_symbs = fo_odesys.dep_var_func_symbs

        # Init attributes for possible analytically solvable y's
        self._init_val_symbs = {}
        self._solved = OrderedDict()


    def mk_init_val_symb(self, y):
        new_symb = sympy.symbols(str(y) + '_init')
        assert not new_symb == self._fo_odesys.indep_var_symb
        assert not new_symb in self._fo_odesys.dep_var_func_symbs
        assert not new_symb in self._fo_odesys.param_symbs
        return new_symb


    def recursive_analytic_reduction(self):
        """
        Attempts to solve some y's analytically

        TODO: recreate possible 2nd order ODE and check solvability
        """
        x = self._fo_odesys.indep_var_symb
        solved = {}
        new_init_val_param_symbs = {}
        changed_last_loop = True
        while changed_last_loop:
            changed_last_loop = False
            for yi, expr in self._fo_odesys.f.iteritems():
                if yi in solved: continue
                # Check to see if it only depends on itself
                expr = expr.subs(solved)
                Jacobian_row_off_diag = [expr.diff(m) for m \
                                         in self._fo_odesys.dep_var_func_symbs if m != yi]
                if all([c == 0 for c in Jacobian_row_off_diag]):
                    # Attempt solution (actually: assume success)
                    rel = sympy.Eq(yi.diff(x), expr)
                    sol = sympy.dsolve(rel, yi)
                    # Assign new symbol to inital value
                    new_symb = self.mk_init_val_symb(yi)
                    new_init_val_param_symbs[yi] = new_symb
                    sol_init_val = determine_const_val_for_init_val(
                        sol, new_symb, x)
                    solved[yi] = sol_init_val.rhs
                    changed_last_loop = True

        # Post processing of solutions
        self._solved.update(solved)
        if len(solved) > 0:
            new_dep_var_func_symbs = [x for x in self._fo_odesys.dep_var_func_symbs \
                                      if x not in solved]
            new_f = {}
            for k in new_dep_var_func_symbs:
                new_f[k] = self._fo_odesys[k].subs(solved)

            self._fo_odesys.dep_var_func_symbs = new_dep_var_func_symbs
            self._fo_odesys.f = new_f

        self._init_val_symbs.update(new_init_val_param_symbs)
        return new_init_val_param_symbs


    def update_init_vals(self, init_vals):
        """
        Updates initial values
        """
        for k, v in initial_values.iteritems:
            if k in self._init_vals:
                self._init_vals[k] = v
            else:
                raise KeyError('Initial value for {} does not exist'.format(k))


    def integrate(self, t0, tend, N, h = None, order = 2):
        """
        Integrates the non-analytic odesystem and evaluates the analytic functions
        for the dependent variables (if there are any).
        """
        if len(self._solved) > 0:
            # If we have solved parts analytically
            a_y0 = {}
            for yi, init_val_symb in self._init_val_symbs.iteritems():
                a_y0[init_val_symb] = self._init_vals.pop(yi)
            self._analytic_evalr = self._AnalyticEvalr(
                self._solved.values(), self._fo_odesys.indep_var_symb,
                self._fo_odesys.param_vals_by_symb, order = order)

        if len(self._init_vals) > 0:
            # If there are any non-analytic equations left
            self._integrator = self._Integrator(self._fo_odesys)
            self._integrator.integrate(self._init_vals, t0, tend, N, h, order)
            self.tout = self._integrator.tout
        else:
            if N == 0: N = self.default_N
            self.tout = np.linspace(t0, tend, N)

        if len(self._solved) > 0:
            self._analytic_evalr.eval_for_indep_array(self.tout, a_y0)

    @property
    def yout(self):
        """
        Unified the output of the numerical and analyitc results.
        """
        nt = self.tout.shape[0]
        _yout = np.empty((nt, len(self._ori_dep_var_func_symbs)), self._dtype)
        j = 0 # Counter of analytic results
        k = 0 # Counter of numerical results
        for i, yi in enumerate(self._ori_dep_var_func_symbs):
            if yi in self._solved:
                _yout[:, i] = self._analytic_evalr.yout[:, j]
                j += 1
            else:
                if len(self._fo_odesys.dep_var_func_symbs) > 0:
                    _yout[:, i] = self._integrator.yout[:, k]
                    k += 1
        return _yout


    def Dy(self):
        nt = self.tout.shape[0]
        _out = np.empty((nt, len(self._ori_dep_var_func_symbs)), self._dtype)
        if len(self._fo_odesys.dep_var_func_symbs) > 0:
            num_dyout = np.array(
                [self._fo_odesys.dydt(
                    t, self._integrator.yout[i,:], self._fo_odesys.params_val_lst) for\
                 (i,), t in np.ndenumerate(self.tout)])
        j = 0 # Counter of analytic results
        k = 0 # Counter of numerical results
        for i, yi in enumerate(self._ori_dep_var_func_symbs):
            if yi in self._solved:
                _out[:, i] = self._analytic_evalr.dyout[:, j]
                j += 1
            else:
                if len(self._fo_odesys.dep_var_func_symbs) > 0:
                    _out[:, i] = num_dyout[:, k]
                    k += 1
        return _out


    def DDy(self):
        nt = self.tout.shape[0]
        _out = np.empty((nt, len(self._ori_dep_var_func_symbs)), self._dtype)
        if len(self._fo_odesys.dep_var_func_symbs) > 0:
            num_ddyout = np.array(
                [self._fo_odesys.d2ydt2(
                    t, self.yout[i,:], self._fo_odesys.params_val_lst) for\
                 (i,), t in np.ndenumerate(self.tout)])
        j = 0 # Counter of analytic results
        k = 0 # Counter of numerical results
        for i, yi in enumerate(self._ori_dep_var_func_symbs):
            if yi in self._solved:
                _out[:, i] = self._analytic_evalr.ddyout[:, j]
                j += 1
            else:
                if len(self._fo_odesys.dep_var_func_symbs) > 0:
                    _out[:, i] = num_ddyout[:, k]
                    k += 1
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
                     in self._ori_dep_var_func_symbs])[str(symb)][:, 0]

    def plot(self, indices = None, interpolate = False, show = True):
        """
        Rudimentary plotting utility for quick inspection of solutions
        TODO: move this from here,  make more general to accept mixed ODE sol +
        analytic y curves
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
            lbl = str(self._ori_dep_var_func_symbs[i])
            if interpolate:
                plt.plot(ipx, ipy[:, i], label = lbl + ' (interpol.)',
                         marker = 'None', ls = lsi, color = ci)
                lsi = 'None'
            plt.plot(self.tout, self.yout[:, i], label = lbl,
                     marker = mi, ls = lsi, color = ci)
            plt.plot()
        plt.legend()
        if show: plt.show()
