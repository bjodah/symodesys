#!/usr/bin/env python
# -*- coding: utf-8 -*-

# other imports
import numpy as np

def get_new_symbs(expr, known_symbs):
    new_symbs = set()
    for atom in expr.atoms():
        if not atom in known_symbs and not atom.is_Number:
            new_symbs.add(atom)
    return new_symbs


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
