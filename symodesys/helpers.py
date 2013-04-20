
# stdlib imports
import os
import imp
from functools import wraps
from collections import OrderedDict # for OrderedDefaultdict
from hashlib import md5

# other imports
import sympy
import numpy as np
from sympy.utilities.autowrap import autowrap, ufuncify

def cache(f):
    data = {}
    @wraps(f)
    def wrapper(*args):
        hashable_args=[]
        for x in args:
            if isinstance(x, dict):
                hashable_args.append(frozenset(x.items()))
            elif isinstance(x, list):
                hashable_args.append(tuple(x))
            else:
                hashable_args.append(x)
        hashable_args = tuple(hashable_args)
        if not hashable_args in data:
            data[hashable_args] = f(*args)
        return data[hashable_args]
    wrapper.cache_clear = lambda: data.clear()
    wrapper.cache = data
    return wrapper

def deprecated(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        print('This is a deprecation warning regarding: {}'.format(
            f.__name__))
        return f(*args, **kwargs)
    return wrapper


def md5_of_file(path):
    md = md5()
    with open(path,'rb') as f:
        for chunk in iter(lambda: f.read(128*md.block_size), b''):
             md.update(chunk)
    return md.digest()

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

    def __init__(self, order=0, dtype=None):
        """

        Arguments:
        - `exprs`: List of expressions to be evaluated
        - `indep_var_symb`: Sympy symbol of the independent variable
        - `params_by_symb`: Dictionary mapping sympy symbols of parameters to values
        - `order`: set higher than 0 (default) in order to also evaluate derivatives.
                   (output is sotred in self.Yout)
        """
        self.order = order
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
        _Yout = np.empty((len(arr), len(self._exprs), self.order+1), dtype=self._dtype)
        for expr_idx, expr in enumerate(self._exprs):
            for order in range(self._order):
                subs = self._params_by_symb
                subs.update(extra_params_by_symb)
                for i, t in enumerate(arr):
                    subs.update({self._indep_var_symb: t})
                    diff_expr = expr.diff(self._indep_var_symb, order)
                    _Yout[i, expr_idx, order] = diff_expr.subs(subs)
        self.Yout = _Yout


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
