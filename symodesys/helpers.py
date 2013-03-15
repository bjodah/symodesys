
# stdlib imports
import imp
from functools import wraps
from collections import OrderedDict # for OrderedDefaultdict

# other imports
import sympy
import numpy as np
from sympy.utilities.autowrap import autowrap

def cache(f):
    data = {}
    @wraps(f)
    def wrapper(*args):
        if not args in data:
            data[args] = f(*args)
        return data[args]
    return wrapper

def subs_set(s, subsd):
    """
    Substititues entities in a set ``s'' for matching keys in ``subsd''
    with corresponding values
    """
    t = s.copy()
    for k, v in subsd.iteritems():
        if k in s:
            s.remove(k)
            s.add(v)
    return t

class OrderedDefaultdict(OrderedDict):
    """
    From http://stackoverflow.com/questions/4126348/\
    how-do-i-rewrite-this-function-to-implement-ordereddict/4127426#4127426
    """

    def __init__(self, *args, **kwargs):
        newdefault = None
        newargs = ()
        if args:
            newdefault = args[0]
            if not (newdefault is None or callable(newdefault)):
                raise TypeError('first argument must be callable or None')
            newargs = args[1:]
        self.default_factory = newdefault
        super(self.__class__, self).__init__(*newargs, **kwargs)

    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):  # optional, for pickle support
        args = self.default_factory if self.default_factory else tuple()
        return type(self), args, None, None, self.items()


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

def import_(filename):
    """ Imports (cython generated) .so file """
    path, name = os.path.split(filename)
    name, ext = os.path.splitext(name)
    fobj, filename, data = imp.find_module(name, [path])
    mod = imp.load_module(name, fobj, filename, data)
    return mod

def get_new_symbs(expr, known_symbs):
    new_symbs = set()
    for atom in expr.atoms():
        if not atom in known_symbs and not atom.is_Number:
            new_symbs.add(atom)
    return new_symbs


class PieceWiseShiftedPolyTraj(object):
    """
    Class to store a series of piece-wise
    shifted polynomials

    Note: this does not do fitting, it only
          solves for one single unique solution
          Hence possible orders include: 1,3,5,(7),..
          were 7 and above needs entries in A
    """

    # The coefficients of the polynomial is determined
    # for shifted and scaled interval (x: {0,1}) data
    # is ordered as {y(0),y(1),y'(0),y'(1),y''(0),y''(1)}
    A = {1: np.array([
            [1, 0],
            [1, 1]], dtype = np.float64),
         3: np.array([
            [1, 0, 0, 0],
            [1, 1, 1, 1],
            [0, 1, 0, 0],
            [0, 1, 2, 3]], dtype = np.float64),
         5: np.array([
            [1, 0, 0, 0,  0,  0],
            [1, 1, 1, 1,  1,  1],
            [0, 1, 0, 0,  0,  0],
            [0, 1, 2, 3,  4,  5],
            [0, 0, 2, 0,  0,  0],
            [0, 0, 2, 6, 12, 20]], dtype = np.float64),
            }

    _xsymb = sympy.Symbol('x', real=True)
    @property
    def _coeffsymb(self):
        [sympy.Symbol('c_' + str(i), real=True) for i in range(self.order + 1)]

    def __init__(self, t, Y):
        """
        t = parameter data points
        Y.shape
        len(Y[i][0]) === order+1
        """
        if isinstance(t, np.ndarray):
            self.t = t
        else:
            self.t = np.array(t)
        if isinstance(Y, np.ndarray):
            self.Y = Y
        else:
            self.Y = np.array(Y)
        assert self.Y.shape[0] == len(self.t) - 1
        self.order = self.Y.shape[1] * 2 - 1
        self._fit_poly()
        self._mk_eval()

    def _fit_poly(self):
        self._coeff_shifted_scaled = np.empty(self.Y.shape)
        b = np.empty((self.order + 1, 1))
        for i in range(len(t)-1):
            for j in range(self.order / 2 + 1):
                b[j * 2] = self.Y[i, j]
                b[j * 2 + 1] = self.Y[i + 1, j]
            self._coeff_shifted_scaled[i,:] = np.linalg.solve(A[self.order], b)

    @property
    def _coeff_shifted(self):
        dt = np.diff(self.t)
        k = 1 / dt
        coeff_factor = np.vander(k, self.order + 1)
        return self._coeff_shifted_scaled * coeff_factor[:,::-1]

    @property
    def _coeff(self):
        conv_cb, conv_args = self._mk_converter()
        coeff = np.empty(self._coeff_shifted.shape)
        for i in range(self._coeff_shifted.shape[1]):
            coeff[:, i] = conv_vb[i](
                self.t[:-1], self._coeff_shifted[:, conv_args[i]])
        return coeff

    @property
    def _base_poly(self):
        return sum([self._coeffsymb[i]*self._xsymb**i for i in range(self.order + 1)])

    @cache
    def _mk_converter(self):
        c = self._coeffsymb
        x = self._xsymb
        x0 = sympy.Symbol('x0', real=True)
        p = self._base_poly.subs({x: x - x0}).expand()
        k = []
        for o in reversed(range(1, self.order + 1)):
            c = p.coeff(x ** o)
            k.append(c)
            p = sum([e for e in p.args if e not in c*x**o])
        k.appned(p)
        k = reversed(k)
        bin_args = [[x0] + [ci for ci in c if ci in ki] for ki in k]
        bin_cb = [ufuncify(bin_args[i], ki) for i, ki in enumerate(k)]
        return bin_cb, [[i for i in len(c) if c[i] in ki] for ki in k]

    def _mk_eval(self):
        coeff = self._coeff
        poly = self._base_poly.subs(dict(zip(self._coeffsymb, coeff)))
        self._eval = ufuncify([self._xsymb] + self._coeffsymb, poly)

    def __call__(self, t):
        return self._eval(t)
