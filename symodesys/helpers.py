from functools import wraps

import numpy as np

def cache(f):
    data = {}
    @wraps(f)
    def wrapper(*args):
        if args in data:
            return data[args]
        else:
            data[args] = f(*args)
            return data[args]
    return wrapper

class SympyEvalr(object):
    """
    Evalurates sympy expressions dependent on one variable
    also substitutes parameters in expression using provided dict
    Also evaluates the derivates in the independent variable up to
    requested order.
    """

    default_dtype = np.float64

    def __init__(self, exprs, indep_var_symb, params_by_symb, order=0,
                 dtype=None):
        """

        Arguments:
        - `exprs`: List of expressions to be evaluated
        - `indep_var_symb`: Sympy symbol of the independent variable
        - `params_by_symb`: Dictionary mapping sympy symbols of parameters to values
        - `order`: set higher than 0 (default) in order to also evaluate derivatives.
                   (output is sotred in self.yout, self.dyout, self.ddyout...  )
        """
        self._exprs = exprs
        self._indep_var_symb = indep_var_symb
        self._params_by_symb = params_by_symb
        self._order = order
        if dtype == None: dtype = self.default_dtype
        self._dtype = dtype

    def _init_out_arrays(self, n):
        self._out = {}
        for i in range(self._order + 1):
            self._out[i] = np.empty((n, len(self._exprs)), dtype=self._dtype)

    def __getattr__(self, attr):
        if attr.endswith('yout'):
            dstr = attr[:-4]
            n = len(dstr)
            if dstr != 'd' * n: raise AttributeError
            return self._out[n]
        else:
            raise AttributeError("%r object has no attribute %r" %
                                 (type(self).__name__, attr))

    def eval_for_indep_array(self, arr, extra_params_by_symb):
        """
        Evaluate all expressions for values of indepedndent variable
        in array `arr` using self._params_symb and `extra_params_by_symb`
        for static substitution in sympy expressions in the list self._exprs
        """
        self._init_out_arrays(len(arr))
        for expr_idx, expr in enumerate(self._exprs):
            for order in range(self._order):
                subs = self._params_by_symb
                subs.update(extra_params_by_symb)
                for i, t in enumerate(arr):
                    subs.update({self._indep_var_symb: t})
                    diff_expr = expr.diff(self._indep_var_symb, order)
                    self._out[order][i, expr_idx] = diff_expr.subs(subs)


def plot_numeric_vs_analytic(Sys, indep_var_lim,
                             init_dep_var_vals_by_token, param_vals, N = 0):
    """
    Integrate
    """
    sys = Sys()
    sys.update_params_by_token(param_vals)

    y0 = {sys[k]: v for k, v in init_dep_var_vals_by_token.items()}
    ivp = IVP(sys, y0)

    t0, tend = indep_var_lim

    ivp.integrate(t0, tend, N)
    t, y = ivp.tout, ivp.yout

    plt.subplot(311)
    ivp.plot(interpolate = True, show = False)

    analytic_y = sys.analytic_y(t, init_dep_var_vals_by_token)
    plt.subplot(312)
    plt.plot(t, (y[:, 0] - analytic_y) / ivp._integrator.abstol,
             label = 'abserr / abstol')
    plt.legend()
    plt.subplot(313)
    plt.plot(t, (y[:, 0] - analytic_y) / analytic_y / ivp._integrator.reltol,
             label = 'relerr / reltol')
    plt.legend()
    plt.show()
