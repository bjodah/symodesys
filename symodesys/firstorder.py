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

    #
    param_tokens = None
    _param_basesymb = 'k'
    param_default_values = None

    #_attrs_to_cmp is used for checking equality of class instances
    _attrs_to_cmp = ['indep_var_symb', 'num_dep_vars', 'num_params',
                     'dep_var_func_symbs', 'param_symbs', 'f']

    def __eq__(self, other):
        for attr in self._attrs_to_cmp:
            if getattr(self, attr) != getattr(other, attr): return False
        return True

    def __init__(self, params_by_tokens = None):
        self._init_dep_var_tokens()
        self._init_param_tokens()
        if params_by_tokens == None:
            params_by_tokens = {}
        else:
            for token in params_by_tokens:
                assert token in self.param_tokens

        if self.param_default_values == None:
            self._param_default_values = dict(
                [(self[k], params_by_tokens.get(k, 0.0)) for k in self.param_tokens])


    def _init_dep_var_tokens(self):
        if self.dep_var_tokens == None:
            self.dep_var_tokens = [self._dep_var_basesymb + str(i) for i in range(self.num_dep_vars)]
        assert len(self.dep_var_tokens) == len(set(self.dep_var_tokens))

    def _init_param_tokens(self):
        if self.param_tokens == None:
            self.param_tokens = [self._param_basesymb + str(i) for i in range(self.num_params)]
        assert len(self.dep_var_tokens) == len(set(self.dep_var_tokens))

    def __getitem__(self, key):
        """
        If one wants to access the symbol of a dep_var_func_symbs or a param_symbs
        and do not want to hardcode the order in the code for item access, it can be retrieved
        using this function
        """
        if self.param_tokens == None: self._init_param_tokens()
        if self.dep_var_tokens == None: self._init_dep_var_tokens()

        if key in self.dep_var_tokens:
            assert key not in self.param_tokens
            return sympy.Function(key)(self.indep_var_symb)
        elif key in self.param_tokens:
            return sympy.symbols(key)
        else:
            raise KeyError('Unknown token')


    @property
    def dep_var_func_symbs(self):
        """
        May be subclassed
        should return list of sympy.Function(``token_string'')(self.indep_var) instances
        The order in this list defines indices in vectors and matrices used by underlying
        numerical integration.
        """
        #return sympy.symarray(self._dep_var_basesymb, self.num_dep_vars)
        # else:
        #     tokens = self.dep_var_tokens
        #     # Make sure there are no duplicates
        #     assert len(tokens) == len(set(tokens))

        #return [sympy.Function(y)(self.indep_var_symb) for y in self.dep_var_tokens]
        return [self[y] for y in self.dep_var_tokens]


    @property
    def param_symbs(self):
        """
        To be subclassed
        should return e.g.  sympy.symarray('k', 7)
        """
        return [self[k] for k in self.param_tokens]



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
            lambda q, i: self.f[self.dep_var_func_symbs[i]]
            )

    def _get_all_num_subs(self, indep_val, dep_vals, param_vals):
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
        all_subs = self._get_all_num_subs(indep_val, dep_vals, param_vals)
        return [x.subs(all_subs) for x in [self.f[k] for k in self.dep_var_func_symbs]]

    def dydt_jac(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating jacobian of dydt for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self.dep_var_func_symbs, self._params_symbs
        for provided values of the independent, the dependent and parameter
        variables (provided as dictionaries)
        """
        all_subs = self._get_all_num_subs(indep_val, dep_vals, param_vals)
        #all_subs = dict(indep_val.items() + dep_vals.items() + param_vals.items()
        return [[cell.subs(all_subs) for cell in row] for row \
                in self._fmat().jacobian(self.dep_var_func_symbs).tolist()]


    def d2ydt2(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating d2ydt2 (f) for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self.dep_var_func_symbs, self._params_symbs
        for provided values of the independent, the dependent and parameter
        variables (provided as dictionaries
        """
        f = dict([(k.diff(self.indep_var_symb), v) for k, v in self.f.iteritems()])
        dfdt_lst = [self.f[y].diff(self.indep_var_symb).subs(f) for y in self.dep_var_func_symbs]
        all_num_subs = self._get_all_num_subs(indep_val, dep_vals, param_vals)
        return [dfdt.subs(all_num_subs) for dfdt in dfdt_lst]


        # partial_f_partial_t = sympy.Matrix(
        #     self.num_dep_vars, 1,
        #     lambda i, q: self.f[self.dep_var_func_symbs[i]].diff(self.indep_var_symb)
        #     )
        # # We only want expressions explicitly dependent on t
        # partial_f_partial_t = partial_f_partial_t.subs(
        #     dict(zip([x.diff(self.indep_var_symb) for x in self.dep_var_func_symbs],
        #              [0] * self.num_dep_vars)))
        # d2ydt2_expr = self._fmat().jacobian(self.dep_var_func_symbs) * \
        #                self._fmat().transpose() + partial_f_partial_t
        # sympy.pprint(d2ydt2_expr)

        # return [x.subs(self._get_all_num_subs(indep_val, dep_vals, param_vals)) for x \
        #         in d2ydt2_expr]

    def transform_indep_var_to_log_scale(self):
        pass

    def transform_dep_vars_to_log_scale(self):
        pass

    def reduce_sys_by_solving_decoupled_vars_analytically(self):
        pass
