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

    indep_var_symb = None # ODE implies 1 indep. variable, set to sympy.symbol(...)
    dep_var_func_symbs = None
    param_symbs = None
    f = None
    param_vals_by_symb = None


    #_attrs_to_cmp is used for checking equality of class instances
    _attrs_to_cmp = ['indep_var_symb', 'dep_var_func_symbs', 'param_symbs',
                     'f', 'param_vals_by_symb']

    def __init__(self):
        # By use of property functions makes overwriting harder
        # therefore 'indep_var_symb', 'dep_var_func_symbs', 'param_symbs', 'f'
        # are initialized via instance methods that by default initializes
        # empty list / dict only if attribute is missing
        self._init_dep_var_func_symbs()
        self._init_param_symbs()
        self.init_f()

    def _init_dep_var_func_symbs(self):
        """
        To be subclassed (or add list prop: dep_var_func_symbs)

        should return list of sympy.Function(``token_string'')(self.indep_var) instances
        The order in this list defines indices in vectors and matrices used by underlying
        numerical integration.
        """
        if self.dep_var_func_symbs == None:
            self.dep_var_func_symbs = []

    def _init_param_symbs(self):
        """
        To be subclassed (or add list prop: param_symbs)

        should return list of sympy.symbols(``token_sting'') instances
        The order in this list defines indices in vectors and matrices used by underlying
        numerical integration.
        (When subclassing, sympy.symarray might be useful.)
        """
        if self.param_symbs == None:
            self.param_symbs = []

    def init_f(self):
        """
        To be subclassed (or add dict prop: f)

        self.f should return a dict of length
        len(self.dep_var_func_symb) for the first-order derivatives
        of the self.dep_var_func_symbs (the underived dep_var_func_symb acts as key)
        expressed solely in numerical constants, sympy function expressions,
        indep_var_symb, dep_var_func_symbs and param_symbs"""
        if self.f == None:
            self.f = {}

    def _init_param_vals_by_symb(self):
        """
        To be subclassed (or add dict prop: param_vals_by_symb)

        should return a dict of length
        len(self.param_symbs) with the numerical values for each

        The values of the parameters are included in this low level
        both for programmatical convenience (not needing to pass values
        through between interfaces), but it can also be viewed as part
        of the ODE system (especially considering special cases when
        some parameter is 0)
        """
        if self.param_vals_by_symb == None:
            self.param_vals_by_symb = {}

    def __eq__(self, other):
        for attr in self._attrs_to_cmp:
            if getattr(self, attr) != getattr(other, attr): return False
        return True

    def _fmat(self):
        """
        Transorms the dict mapping dep_var_func_symbs self.f
        """
        return sympy.Matrix(
            1, len(self.dep_var_func_symbs),
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

    @property
    def params_val_lst(self):
        return [self.param_vals_by_symb[k] for k in self.param_symbs]


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

    indep_var_symb = sympy.symbols('t') # ODE implies 1 indep. variable

    # The following must be provided (string reprs instead of sympy.symbols)
    dep_var_tokens = None
    param_tokens = None


    params_by_token = None


    def __init__(self):
        if self.param_tokens == None:
            self.param_tokens = []
        if self.params_by_token == None:
            self.params_by_token = {}
        super(SimpleFirstOrderODESystem, self).__init__()

    def update_params_by_token(self, params_by_token):
        for token in params_by_token:
            if not token in self.param_tokens: raise KeyError(
                'Parameter token ``{}" unknown'.format(token))

        if self.params_by_token == None:
            self.params_by_token = {}
        self.params_by_token.update(params_by_token)

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
    def param_vals_by_symb(self):
        return dict([(self[k], self.params_by_token[k]) for\
                     k in self.param_tokens])

    def _init_dep_var_func_symbs(self):
        assert self.dep_var_func_symbs == None
        self.dep_var_func_symbs = [self[y] for y in self.dep_var_tokens]

    def _init_param_symbs(self):
        assert self.param_symbs == None
        self.param_symbs = [self[k] for k in self.param_tokens]

    # End overloading
