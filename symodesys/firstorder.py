from symodesys.odesys import ODESystem

from operator import or_
from collections import OrderedDict

import sympy


class FirstOrderODESystem(ODESystem):
    """
    When the ODE systems equations is generated from user
    data this class is to be subclassed to provide the routines
    described below
    """
    # TODO: implement the routines for variable substitution

    #_attrs_to_cmp_for_eq is used for checking equality of class instances
    _attrs_to_cmp_for_eq = ['indep_var_symb', 'dep_var_func_symbs',
                            'param_symbs', 'f']

    def __init__(self):
        # Using property functions makes overwriting harder therefore
        # 'indep_var_symb', 'dep_var_func_symbs', 'param_symbs', 'f'
        # are initialized via instance methods that by default initializes
        # empty list / dict only if attribute is missing
        self._init_dep_var_func_symbs()
        self._init_param_symbs()
        self.init_f()

    @property
    def _odeqs_by_indep_var(self):
        return OrderedDict((k, (1, self.f[k])) for k \
                           in self.dep_var_func_symbs)

    def _init_dep_var_func_symbs(self):
        """
        To be subclassed (or add list prop: dep_var_func_symbs)

        should return list of
          sympy.Function(``token_string'')(self.indep_var) instances
        The order in this list defines indices in vectors and matrices
        used by underlying numerical integration.
        """
        if self.dep_var_func_symbs == None:
            self.dep_var_func_symbs = []

    def _init_param_symbs(self):
        """
        To be subclassed (or add list prop: param_symbs)

        should return list of sympy.symbols(``token_sting'') instances
        The order in this list defines indices in vectors and matrices
        used by underlying numerical integration.
        (When subclassing, sympy.symarray might be useful.)
        """
        if self.param_symbs == None:
            self.param_symbs = []

    def init_f(self):
        """
        To be subclassed (or add dict prop: f)

        self.f should return a dict of length
        len(self.dep_var_func_symb) for the first-order derivatives
        of the self.dep_var_func_symbs (the underived dep_var_func_symb
        acts as key) expressed solely in numerical constants, sympy
        function expressions, indep_var_symb, dep_var_func_symbs and
        param_symbs"""
        if self.f == None:
            self.f = {}

    # NOTE: we don't take ODESystem to be immutable anymore,
    #       maybe we need an hashable export?
    # def __hash__(self):
    #     """ ODESystem is taken to be immutable """
    #     hashes = [hash(x) for x in self._attrs_to_cmp if x != 'f']
    #     fhash = hash(frozenset(self.f.items()))
    #     return sum(hashes) + fhash

    def _fmat(self):
        """
        Transforms the dict mapping dep_var_func_symbs self.f
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

    def param_val_lst(self, param_vals_by_symb):
        return [param_vals_by_symb[k] for k in self.param_symbs]


    def dydt(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating dydt (f) for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self.dep_var_func_symbs, self._param_symbs
        for provided values of the independent, the dependent and
        parameter variables (with same order)

        Note: The signature of the function employs float point data
              (or lists thereof) in order to be compatible with
              e.g. scipy integrators, hence _get_all_numb_subs
        """
        all_subs = self._get_all_num_subs(indep_val, dep_vals, param_vals)
        return [x.subs(all_subs) for x in \
                [self.f[k] for k in self.dep_var_func_symbs]]


    def dydt_jac(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating jacobian of dydt for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self.dep_var_func_symbs,
        self._param_symbs for provided values of the independent,
        the dependent and parameter
        variables (provided as dictionaries)

        Note: The signature of the function employs float point
               data (or lists thereof) in order to be compatible
               with e.g. scipy integrators, hence _get_all_numb_subs
        """
        all_subs = self._get_all_num_subs(
            indep_val, dep_vals, param_vals)
        return [[cell.subs(all_subs) for cell in row] for row in \
                self._fmat().jacobian(self.dep_var_func_symbs).tolist()]


    def d2ydt2(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating d2ydt2 (f) for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self.dep_var_func_symbs,
        self._param_symbs for provided values of the independent,
        the dependent and parameter variables (provided as dictionaries

        Note: The signature of the function employs float point
        data (or lists thereof) in order to be compatible with
        e.g. scipy integrators, hence _get_all_numb_subs
        """
        f = {k.diff(self.indep_var_symb): v for k, v in \
             self.f.iteritems()}
        dfdt_lst = [self.f[y].diff(self.indep_var_symb).subs(f) for \
                    y in self.dep_var_func_symbs]
        all_num_subs = self._get_all_num_subs(
            indep_val, dep_vals, param_vals)
        return [dfdt.subs(all_num_subs) for dfdt in dfdt_lst]


    def transform_indep_var_to_log_scale(self):
        # TODO should be more general than just log_scale: variable subst
        pass

    def transform_dep_vars_to_log_scale(self):
        # TODO should be more general than just log_scale: variable subst
        pass


class SimpleFirstOrderODESystem(FirstOrderODESystem):
    """
    This class provides convenience methods for generating the
    symbols of the idependent variable symbol, dependent variable symbols
    and parameter symbols. It is useful when the equations are not
    algorithmatically generated but by user subclassing (of this class).
    """

    indep_var_symb = sympy.symbols('t') # ODE implies 1 indep. variable

    # The following must be provided in this "Simple" subclass
    # (string reprs instead of sympy.symbols)
    dep_var_tokens = None
    param_tokens = None

    def __init__(self):
        if self.param_tokens == None:
            self.param_tokens = []
        super(SimpleFirstOrderODESystem, self).__init__()

    def get_param_vals_by_symb_from_by_token(self, param_vals_by_token):
        for token in param_vals_by_token:
            if not token in self.param_tokens: raise KeyError(
                'Parameter token ``{}" unknown'.format(token))
        return {self[k]: v for k, v in param_vals_by_token.iteritems()}

    def __getitem__(self, key):
        """
        If one wants to access the symbol of a dep_var_func_symbs
        or a param_symbs and do not want to hardcode the order in
        the code for item access, it can be retrieved using this function
        """
        if key in self.dep_var_tokens:
            assert key not in self.param_tokens
            return sympy.Function(key)(self.indep_var_symb)
        elif key in self.param_tokens:
            return sympy.symbols(key)
        else:
            raise KeyError('Unknown token')


    # Begin overloading

    def param_vals_by_symb(self, params_by_token):
        return dict([(self[k], params_by_token[k]) for\
                     k in self.param_tokens])

    def _init_dep_var_func_symbs(self):
        # The assert is there to signal need to subclass if using
        # SimpleFirstOrderODESystem
        assert self.dep_var_func_symbs == None
        self.dep_var_func_symbs = [self[y] for y in self.dep_var_tokens]

    def _init_param_symbs(self):
        # The assert is there to signal need to subclass if using
        # SimpleFirstOrderODESystem
        assert self.param_symbs == None
        self.param_symbs = [self[k] for k in self.param_tokens]

    # End overloading
