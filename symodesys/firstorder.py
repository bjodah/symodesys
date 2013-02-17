import sympy

from operator import or_

class FirstOrderODESystem(object):
    """
    When the ODE systems equations is generated from user
    data this class is to be subclassed to provide the routines
    described below

    TODO: implement the routines for variable substitution
    TODO: add properties(?) for is_autonomous and is_linear
    """

    indep_var_symb = sympy.symbols('t') # ODE implies 1 indep. variable

    num_dep_vars = None # Size of ODE system
    num_params = None # Number of parameters

    #_attrs_to_cmp is used for checking equality of class instances
    _attrs_to_cmp = ['indep_var_symb', 'dep_var_func_symbs', 'param_symbs', 'f']

    def __eq__(self, other):
        for attr in self._attrs_to_cmp:
            if getattr(self, attr) != getattr(other, attr): return False
        return True

    @property
    def default_params(self):
        pass


    @property
    def dep_var_func_symbs(self):
        """
        To be subclassed.
        should return list of sympy.Function(``token_string'')(self.indep_var) instances
        The order in this list defines indices in vectors and matrices used by underlying
        numerical integration.
        """
        pass

    @property
    def param_symbs(self):
        """
        To be subclassed.
        should return list of sympy.symbols(``token_sting'') instances
        The order in this list defines indices in vectors and matrices used by underlying
        numerical integration.
        (When subclassing, sympy.symarray might be useful.)
        """
        pass

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

        Note: The signature of the function employs float point data (or lists thereof)
        in order to be compatible with e.g. scipy integrators, hence _get_all_numb_subs
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

        Note: The signature of the function employs float point data (or lists thereof)
        in order to be compatible with e.g. scipy integrators, hence _get_all_numb_subs
        """
        all_subs = self._get_all_num_subs(indep_val, dep_vals, param_vals)
        return [[cell.subs(all_subs) for cell in row] for row \
                in self._fmat().jacobian(self.dep_var_func_symbs).tolist()]


    def d2ydt2(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating d2ydt2 (f) for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self.dep_var_func_symbs, self._params_symbs
        for provided values of the independent, the dependent and parameter
        variables (provided as dictionaries

        Note: The signature of the function employs float point data (or lists thereof)
        in order to be compatible with e.g. scipy integrators, hence _get_all_numb_subs
        """
        f = dict([(k.diff(self.indep_var_symb), v) for k, v in self.f.iteritems()])
        dfdt_lst = [self.f[y].diff(self.indep_var_symb).subs(f) for \
                    y in self.dep_var_func_symbs]
        all_num_subs = self._get_all_num_subs(indep_val, dep_vals, param_vals)
        return [dfdt.subs(all_num_subs) for dfdt in dfdt_lst]


    def transform_indep_var_to_log_scale(self):
        # TODO this should be more general than just log_scale: variable subst!
        pass

    def transform_dep_vars_to_log_scale(self):
        # TODO this should be more general than just log_scale: variable subst!
        pass

    @property
    def eqs(self):
        """
        Returns a list of Sympy Eq instances describing the ODE system
        """
        return [sympy.Eq(k.diff(self.indep_var_symb), v) for k, v in self.f.iteritems()]


class SimpleFirstOrderODESystem(FirstOrderODESystem):
    """
    This class provides convenience methods for generating the
    symbols of the idependent variable symbol, dependent variable symbols
    and parameter symbols. It is useful when the equations are not
    algorithmatically generated but by user subclassing (of this class).
    """

    # if param_symbs is overlaoded (or param_tokens provided)
    # num_params is set in self._init_param_tokens, else it must be set here:

    # Overload dep_var_tokens to get other names than y0, y1, y2...
    # set to e.g. ['f', 'g', 'h']
    # This is useful when manually initializing a FirstOrderODESystem
    # instance, for algorithmic generation, tokens need not be used
    # (overload symbs attributes instead)
    dep_var_tokens = None
    _dep_var_basesymb = 'y'

    #
    param_tokens = None
    _param_basetoken = 'k'
    default_params_by_token = None


    def __init__(self):
        self._init_dep_var_tokens()
        self._init_param_tokens()

    def update_default_params_by_token(self, params_by_token):
        for token in params_by_token:
            if not token in self.param_tokens: raise KeyError(
                'Parameter token ``{}" unknown'.format(token))

        if self.default_params_by_token == None:
            self.default_params_by_token = {}
        self.default_params_by_token.update(params_by_token)

    def _init_dep_var_tokens(self):
        if self.dep_var_tokens == None:
            self.dep_var_tokens = [self._dep_var_basesymb + str(i) for\
                                   i in range(self.num_dep_vars)]
        else:
            if self.num_dep_vars == None:
                self.num_dep_vars = len(self.dep_var_tokens)
            else:
                assert self.num_dep_vars == len(self.dep_var_tokens)
        assert len(self.dep_var_tokens) == len(set(self.dep_var_tokens))

    def _init_param_tokens(self):
        if self.param_tokens == None:
            self.param_tokens = [self._param_basesymb + str(i) for\
                                 i in range(self.num_params)]
        else:
            if self.num_params == None:
                self.num_params = len(self.param_symbs)
            else:
                assert self.num_params == len(self.param_symbs)
        assert len(self.dep_var_tokens) == len(set(self.dep_var_tokens))

    def __getitem__(self, key):
        """
        If one wants to access the symbol of a dep_var_func_symbs or a param_symbs
        and do not want to hardcode the order in the code for item access, it can
        be retrieved using this function
        """

        if key in self.dep_var_tokens:
            assert key not in self.param_tokens
            return sympy.Function(key)(self.indep_var_symb)
        elif key in self.param_tokens:
            return sympy.symbols(key)
        else:
            raise KeyError('Unknown token')



    # Begin overloading
    @property
    def default_params_by_token(self):
        return dict([(self[k], self.default_params_by_token[k]) for\
                     k in self.param_tokens])

    @property
    def dep_var_func_symbs(self):
        """
        May be subclassed
        should return list of sympy.Function(``token_string'')(self.indep_var) instances
        The order in this list defines indices in vectors and matrices used by underlying
        numerical integration.
        """
        return [self[y] for y in self.dep_var_tokens]


    @property
    def param_symbs(self):
        """
        May be subclassed
        should return list of sympy.symbols(``token_sting'') instances
        The order in this list defines indices in vectors and matrices used by underlying
        numerical integration.
        (When subclassing, sympy.symarray might be useful.)
        """
        return [self[k] for k in self.param_tokens]

    # End overloading


