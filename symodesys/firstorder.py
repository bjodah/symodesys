import sympy

from operator import or_

class FirstOrderODESystem(object):
    """
    When the ODE systems equations is generated from user
    data this class is to be subclassed to provide the routines
    described below
    """

    # TODO add properties(?) for is_autonomous and is_linear

    indep_var_symb = sympy.symbols('t') # ODE implies 1 indep. variable

    num_dep_vars = None # Size of ODE system
    num_params = None # Number of parameters

    # Overload dep_var_tokens to get other names than y0, y1, y2...
    # set to e.g. ['f', 'g', 'h']
    dep_var_tokens = None
    _dep_var_basesymb = 'y'
    _param_basesymb = 'k'

    #_attrs_to_cmp is used for checking equality of class instances
    _attrs_to_cmp = ['indep_var_symb', 'num_dep_vars', 'num_params',
                     'dep_var_func_symbs', 'param_symbs', 'f']

    def __eq__(self, other):
        for attr in self._attrs_to_cmp:
            if getattr(self, attr) != getattr(other, attr): return False
        return True

    @property
    def dep_var_func_symbs(self):
        """
        May be subclassed
        should return list of sympy.Function(``token_string'')(self.indep_var) instances
        The order in this list defines indices in vectors and matrices used by underlying
        numerical integration.
        """
        #return sympy.symarray(self._dep_var_basesymb, self.num_dep_vars)
        if self.dep_var_tokens == None:
            tokens = [self._dep_var_basesymb + str(i) for i in range(self.num_dep_vars)]
        else:
            tokens = self.dep_var_tokens
            # Make sure there are no duplicates
            assert len(tokens) == len(set(tokens))
        return [sympy.Function(y)(self.indep_var_symb) for y in tokens]


    @property
    def param_symbs(self):
        """
        To be subclassed
        should return e.g.  sympy.symarray('k', 7)
        """
        return sympy.symarray(self._param_basesymb, self.num_params)

    @property
    def param_default_values(self):
        """
        To be subclassed
        should return list of same length as self.param_symbs()
        with default values. None represents no default value.
        """
        return np.zeros(self.num_params)


    @property
    def f(self):
        """
        To be subclassed, should return a dict of length
        len(self.dep_var_func_symb) for the first-order derivatives
        of the self.dep_var_func_symbs (the underived dep_var_func_symb acts as key)
        expressed solely in numerical constants, sympy function expressions,
        indep_var_symb, dep_var_func_symbs and param_symbs
        """
        pass

    def _fmat(self):
        """
        Transorms the dict mapping dep_var_func_symbs self.f
        """
        return sympy.Matrix(
            1, self.num_dep_vars,
            lambda q, i: self.f[self.dep_var_func_sumbs[i]]
            )

    def _get_all_subs(self, indep_val, dep_vals, param_vals):
        indep_subs = {self.indep_var_symb: indep_val}
        dep_subs = dict(zip(self.dep_var_func_symbs, dep_vals))
        param_subs = dict(zip(self.param_symbs, param_vals))
        all_subs = {}
        for d in [indep_subs, dep_subs, param_subs]:
            all_subs.update(d)
        return all_subs

    def dydt(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating dydt (f) for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self.dep_var_func_symbs, self._params_symbs
        for provided values of the independent, the dependent and parameter
        variables (with same order)
        """
        return [x.subs(self._get_all_subs(indep_val, dep_vals, param_vals)) for x \
                in self.f]

    def dydt_jac(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating jacobian of dydt for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self.dep_var_func_symbs, self._params_symbs
        for provided values of the independent, the dependent and parameter
        variables (with same order)
        """
        all_subs = self._get_all_subs(indep_val, dep_vals, param_vals)
        return [[cell.subs(all_subs) for cell in row] for row \
                in self._fmat().jacobian(self.dep_var_func_symbs).tolist()]


    def d2ydt2(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating d2ydt2 (f) for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self.dep_var_func_symbs, self._params_symbs
        for provided values of the independent, the dependent and parameter
        variables (with same order)
        """
        partial_f_partial_t = sympy.Matrix(
            self.num_dep_vars, 1,
            lambda i, q: self.f[i].diff(self.indep_var_symb)
            )
        d2dydt2_expr = self._fmat().jacobian(self.dep_var_func_symbs) * \
                       self._fmat().transpose() + partial_f_partial_t
        return [x.subs(self._get_all_subs(indep_val, dep_vals, param_vals)) for x \
                in d2dydt2_expr]

    def transform_indep_var_to_log_scale(self):
        pass

    def transform_dep_vars_to_log_scale(self):
        pass

    def reduce_sys_by_solving_decoupled_vars_analytically(self):
        pass
