from symodesys.helpers import OrderedDefaultdict

import sympy

from functools import reduce
from operator import add
from collections import namedtuple

ODEExpr = namedtuple('ODEExpr', 'order expr')


class ODESystem(object):
    """
    Central to the ODESystem if the attribute _odeqs_by_dep_var
    which should be an OrderedDefaultdict(ODEExpr). The keys in the
    dictionary should be of the form:
       sympy.Function('name')(self.indep_var_symb)
    """

    # For tracking what dep var func symbs are generated upon
    # order reduction (instantiate as list)
    # it is useful if one does not want to e.g. plot them
    _1st_ordr_red_helper_fncs = None

    # We need holder for analytically solved parts (instantiate to dict):
    _solved = None

    #_attrs_to_cmp_for_eq is used for checking equality of instances
    _attrs_to_cmp_for_eq = ['indep_var_symb', 'param_symbs']

    indep_var_symb = None # ODE implies 1 indep. variable,
                          #  set to sympy.symbol(...)
    param_symbs = None
    f = None

    dep_var_func_symbs = None

    def __init__(self):
        self._solved = self._solved or {}

    @property
    def known_symbs(self):
        return [self.indep_var_symb] + self._odeqs_by_dep_var.keys() +\
               self.param_symbs + reduce(add, [
            sol_symbs for expr, sol_symbs in self._solved.values()])

    def subs(self, subsd):
        assert all([key in self.known_symbs for key in subsd.keys()])
        self.indep_var_symb = self.indep_var_symb.subs(subsd)
        for key in self._odeqs_by_dep_var:
            self._odeqs_by_dep_var[key].expr = self._odeqs_by_dep_var[key].expr.subs(subsd)
        for i, param in enumerate(self.param_symbs):
            self.param_symbs[i] = self.param_symbs[i].subs(subsd)
        for key, (expr, sol_symbs) in self._solved.iteritems():
            self._solved[key] = expr.subs(subsd), sol_symbs.subs(subsd)

    @property
    def is_first_order(self):
        return all([v.order == 1 for v in self._odeqs_by_dep_var.values()])

    @property
    def is_autonomous(self):
        for dep_var_func_symb, (order, expr) in \
                self._odeqs_by_dep_var.iteritems():
            unfunc_subs = {dvfs: sympy.symbol(dvfs.func.__name__) for \
                           dvfs in self.dep_var_func_symbs}
            if expr.subs(unfunc_subs).diff(self.indep_var_symb) != 0:
                return False
        return True

    @property
    def is_linear(self):
        # Most easily done for first order system?
        for dep_var_func_symb, (order, expr) in \
                self._odeqs_by_dep_var.iteritems():
            for wrt in self.dep_var_func_symbs:
                expr = expr.diff(wrt)

            unfunc_subs = {dvfs: sympy.symbol(dvfs.func.__name__) for \
                           dvfs in self.dep_var_func_symbs}
            if expr.subs(unfunc_subs).diff(self.indep_var_symb) != 0:
                return False
        return True

    @property
    def is_homogeneous(self):
        pass

    def do_sanity_check_of_odeqs(self):
        for fnc, (order, expr) in self._odeqs_by_dep_var.iteritems():
            # Sanity check (indentation level excluded ;-)
            for check_fnc in self._odeqs_by_dep_var.keys():
                # Make sure expr don't contain derivatives larger
                # than in odeqs_by_dep_var order:
                if expr.has(check_fnc):
                    for arg in expr.args:
                        if arg.has(check_fnc):
                            if arg.is_Derivative:
                                fnc, wrt = args[0], args[1:]
                                assert not len(wrt) > \
                                       odeqs_by_dep_var[check_fnc][0]

    def __eq__(self, other):
        for attr in self._attrs_to_cmp_for_eq:
            if getattr(self, attr) != getattr(other, attr): return False
        return True

    @property
    def eqs(self):
        """
        Returns a list of Sympy Eq instances describing the ODE system
        """
        return [self.eq(depvar) for depvar in self._odeqs_by_dep_var.keys()]

    def eq(self, depvar):
        """
        Returns a sympy.Eq for diff eq of ``depvar''
        """
        order, expr = self._odeqs_by_dep_var[depvar]
        return sympy.Eq(depvar.diff(self.indep_var_symb, order), expr)

    @classmethod
    def from_list_of_eqs(cls, lst):
        """
        Determine independnt variable, dependent variables
        and parameter variables and perform Sanity checks
        """
        fncs = []

        # independet variable as key, (order, expr) as value
        odeqs_by_dep_var = OrderedDefaultdict(ODEExpr)
        indep_var_symb = None
        for eq in lst:
            assert eq.lhs.is_Derivative
            fnc, wrt = eq.lhs.args[0], eq.lhs.args[1:]
            assert fnc not in odeqs_by_dep_var
            if indep_var_symb == None:
                assert all([wrt[0] == x for x in wrt]) # No PDEs!
                indep_var_symb = wrt[0]
            else:
                assert all([indep_var_symb == x for x in wrt])
            odeqs_by_dep_var[fnc] = ODEExpr(len(wrt), eq.rhs)

        param_symbs = set()
        for order, expr in odeqs_by_dep_var.values():
            param_symbs.add(get_new_params(expr, odeqs_by_dep_var.keys() + \
                                           [indep_var_symb]))
        new_instance = cls(odeqs_by_dep_var, indep_var_symb,
                           list(param_symbs))
        new_instance.do_sanity_check_of_odeqs()
        return new_instance


def get_new_symbs(expr, known_symbs):
    new_symbs = set()
    for atom in expr.atoms():
        if not atom in known_symbs and not atom.is_Number:
            new_symbs.add(atom)
    return new_symbs


class AnyOrderODESystem(ODESystem):

    _attrs_to_cmp_for_eq = ODESystem._attrs_to_cmp_for_eq +\
                           ['_odeqs_by_dep_var']

    @staticmethod
    def mk_odeqs_by_dep_var():
        """ Convenience function for instantiating OrderedDefaultdict(ODEExpr) """
        return OrderedDefaultdict(ODEExpr)


    # Difference from FirstOrderODESystem
    @property
    def dep_var_func_symbs(self):
        return self._odeqs_by_dep_var.keys()


    def __init__(self, odeqs_by_dep_var, indep_var_symb, param_symbs):
        super(AnyOrderODESystem, self).__init__()
        self._odeqs_by_dep_var = odeqs_by_dep_var
        self.indep_var_symb = indep_var_symb
        self.param_symbs = param_symbs
        self._1st_ordr_red_helper_fncs = self._1st_ordr_red_helper_fncs or []


    def attempt_analytic_sol(self, depvar, hypoexpr, sol_consts):
        """
        Checks if provided analytic (similiar to sympy.solvers.ode.checkodesol)
        expr ``hypoexpr'' solves diff eq of ``depvar'' in odesys. If it does
        new symbols
        """
        eq = self.eq(depvar).subs({depvar: hypoexpr})
        if bool(eq.doit()):
            hypoparams = get_new_symbs(hypoexpr, self.known_symbs)
            self._solved[depvar] = hypoexpr, sol_consts
            return True
        else:
            return False


    @property
    def f(self):
        assert self.is_first_order
        return OrderedDefaultdict(ODEExpr, [\
            (k, v.expr) for k, v in self._odeqs_by_dep_var.items()])

    def get_helper_fnc(self, fnc, order):
        """
        Returns a list of of length order - 1
        for use in reformulation of higher order
        ODE eq in first order ODE (sub) sys.
        """
        helpers = {}
        for i in range(1, order):
            candidate = sympy.Function(str(fnc.func.__name__) + '_h' + \
                                       str(i))(self.indep_var_symb)
            while candidate in self._odeqs_by_dep_var:
                candidate = candidate + '_h'
            helpers[i] = candidate
        return helpers

    def reduce_to_sys_of_first_order(self):
        """ Returns a new instance with reduced order (1st) """
        new_odeqs = self.mk_odeqs_by_dep_var()
        _1st_ordr_red_helper_fncs = self._1st_ordr_red_helper_fncs[:]
        # TODO, revise _1st_ordr_red_helper_fncs format (cur.  len 3 tuple)
        # and update analytic_harmonic_oscillator.py and IVP.plot()
        for fnc, (order, expr) in self._odeqs_by_dep_var.iteritems():
            if order == 1:
                new_odeqs[fnc] = ODEExpr(order, expr)
                continue
            hlpr = self.get_helper_fnc(fnc, order)
            new_odeqs[fnc] = ODEExpr(1, hlpr[1])
            for o in range(1, order - 1):
               new_odeqs[hlpr[o]] = ODEExpr(1, hlpr[o + 1])
            subsd = {sympy.Derivative(fnc, self.indep_var_symb, i): \
                     hlpr[i] for i in range(1, order)}
            new_odeqs[hlpr[order - 1]] = ODEExpr(1, expr.subs(subsd))
            _1st_ordr_red_helper_fncs.extend([
                (fnc, o, expr) for o, expr in hlpr.iteritems()])
        new_instance = self.__class__(new_odeqs, self.indep_var_symb,
                                      self.param_symbs[:])
        new_instance._1st_ordr_red_helper_fncs = _1st_ordr_red_helper_fncs
        return new_instance


class FirstOrderODESystem(ODESystem):
    """
    Special case of ODESystem which the ODESystem is primarily
    given by the attribute `f` which is a dictionary mapping the
    symbol of the dependent variable to the expression of the first
    derivative of same dependent variable with respect to the
    independnt variable.

    When the ODE systems equations is generated from user
    data this class is to be subclassed to provide the routines
    described below
    """
    _attrs_to_cmp_for_eq = ODESystem._attrs_to_cmp_for_eq +\
                           ['f', 'dep_var_func_symbs']

    is_first_order = True # no need to check that in this case

    # TODO: implement the routines for variable substitution

    def __init__(self, odesys = None):
        # Using property functions makes overwriting harder therefore
        # 'indep_var_symb', 'dep_var_func_symbs', 'param_symbs', 'f'
        # are initialized via instance methods that by default initializes
        # empty list / dict only if attribute is missing
        super(AnyOrderODESystem, self).__init__()
        if odesys != None:
            assert odesys.is_first_order
            copy_attrs = ['indep_var_symb', '_1st_ordr_red_helper_fncs',
                          'dep_var_func_symbs', 'param_symbs', 'f', '_solved']
            for attr in copy_attrs:
                setattr(self, attr, getattr(odesys, attr))
        else:
            self._init_dep_var_func_symbs()
            self._init_param_symbs()
            self.init_f()

    @property
    def _odeqs_by_dep_var(self):
        return OrderedDefaultdict(ODEExpr, [(k, ODEExpr(1, self.f[k])) for k \
                           in self.dep_var_func_symbs])

    def recursive_analytic_auto_sol(self):
        """
        Solves equations one by one
        """
        changed_last_loop = True
        # new_y = []
        while changed_last_loop:
            changed_last_loop = False
            for yi, expr in self.f.iteritems():
                if yi in self._solved: continue
                # Check to see if it only depends on itself
                expr = expr.subs(self._solved)
                Jac_row_off_diag = [expr.diff(m) for m \
                                         in self.dep_var_func_symbs if m != yi]
                if all([c == 0 for c in Jac_row_off_diag]):
                    # Attempt solution (actually: assume success)
                    rel = sympy.Eq(yi.diff(self.indep_var_symb), expr)
                    sol = sympy.dsolve(rel, yi)
                    # Assign new symbol to inital value
                    sol_symbs = get_new_symbs(sol.rhs, self.known_symbs)
                    self._solved[yi] = sol_init_val.rhs, sol_symbs
                    changed_last_loop = True
                    # new_y.append(yi)
        # return new_y

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

    Essentially, ``tokens'' are then the generating strings of the
    symbol instances of parent class.
    """

    indep_var_symb = sympy.symbols('t') # ODE implies 1 indep. variable

    # The following must be provided in this "Simple" subclass
    # (string reprs instead of sympy.symbols)
    dep_var_tokens = None
    param_tokens = None

    def __init__(self, *args, **kwargs):
        super(SimpleFirstOrderODESystem, self).__init__(*args, **kwargs)
        self.param_tokens = self.param_tokens or []

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
