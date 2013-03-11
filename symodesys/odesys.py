from symodesys.helpers import OrderedDefaultdict, subs_set, get_new_symbs

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
    _attrs_to_cmp_for_eq = ['indep_var_symb', 'param_symbs', '_solved']

    indep_var_symb = None # ODE implies 1 indep. variable,
                          #  set to sympy.symbol(...)
    param_symbs = None
    f = None


    def __init__(self):
        self._solved = self._solved or {}
        self._1st_ordr_red_helper_fncs = self._1st_ordr_red_helper_fncs or []


    @property
    def non_analytic_depvar(self):
        return [x for x in self._odeqs_by_dep_var if not x in self._solved]

    @property
    def analytic_depvar(self):
        return [x for x in self._odeqs_by_dep_var if x in self._solved]

    @property
    def all_depvar(self):
        return self._odeqs_by_dep_var.keys()

    @property
    def solved_exprs(self):
        """ Convenience attribute """
        return [self._solved[yi][0] for yi in self.analytic_depvar]

    @property
    def known_symbs(self):
        """ Convenience attribute """
        analytic_sol_symbs = set()
        if len(self._solved) > 0:
            for expr, sol_symbs in self._solved.values():
                analytic_sol_symbs = analytic_sol_symbs.union(sol_symbs)
        analytic_sol_symbs = list(analytic_sol_symbs)
        return [self.indep_var_symb] + self.all_depvar +\
               self.param_symbs + analytic_sol_symbs

    def subs(self, subsd):
        # TODO: fix special cases self.indep_var_symb and self.dep_var_func_symbs
        for key in subsd.keys():
            if not key in self.known_symbs:
                raise KeyError('Symbol: {} not known in odesys'.format(key))
        self.indep_var_symb = self.indep_var_symb.subs(subsd)
        for key, odeexpr in self._odeqs_by_dep_var.iteritems():
            self._odeqs_by_dep_var[key] = ODEExpr(
                odeexpr.order, odeexpr.expr.subs(subsd))
        for i, param in enumerate(self.param_symbs):
            self.param_symbs[i] = self.param_symbs[i].subs(subsd)
        for key, (expr, sol_symbs) in self._solved.iteritems():
            self._solved[key] = expr.subs(subsd), subs_set(sol_symbs, subsd)


    def transform_depvars(self, trnsf, inv_subs):
        """
        trnsf: dict mapping old_depvar to tuple of (new_depvar, expr_in_old)
        inv_subs: dict mapping old_depvar to expression in new_depvar
        """
        new_odeqs = OrderedDefaultdict(ODEEq)
        for old_depvar, (order, old_expr) in self._odeqs_by_dep_var.iteritems():
            if not old_depvar in trnsf:
                new_odeqs[old_depvar] = self._odeqs_by_dep_var[old_depvar].subs(inv_subs)
                continue
            # Unpack forward transform
            new_depvar, trnsf_expr = trnsf[old_depvar]
            eqsd = {eq.lhs: eq.rhs for eq in self.eqs}
            # Get new diff eq
            diff_expr = trnsf_expr.diff(self.indep_var_symb, order).subs(eqsd)
            # Save new diff eq (expressed in new depvar through inv_subs)
            new_odeqs[new_depvar] = ODEEq(order, diff_expr.subs(inv_subs))

        # Handle solved:
        new_solved = OrderedDict()
        for old_depvar, old_expr in self._solved.iteritems():
            if not old_depvar in trnsf:
                # Save old analytic expr (expressed in new depvar through inv_subs)
                new_solved[old_depvar] = old_expr.subs(inv_subs)
                continue
            new_depvar, trnsf_expr = trnsf[old_depvar]
            # Save new analytic expr (expressed in new depvar through inv_subs)
            new_solved[new_depvar] = old_expr.subs(inv_subs)
        # Return new instance based on carbon copy of self
        return self.__class__(self, odeqs=new_odeqs, solved=new_solved)

    def transform_indep(self, new_indep_symb, expr_in_old_indep):
        new_odeqs = OrderedDefaultdict(ODEEq)
        for depvar, (order, old_expr) in self._odeqs_by_dep_var.iteritems():
            new_expr = old_expr / sympy.diff(expr_in_old_indep, self.indep_var_symb, order)
            new_expr = new_expr.subs({self.indep_var_symb: sympy.solve(
                new_expr - new_indep_symb, self.indep_var_symb)[0]})
            new_expr = new_expr.subs({expr_in_old_indep: new_indep_symb})
            new_odeqs[depvar] = ODEEq(order, new_expr)

        # Handle solved:
        new_solved = OrderedDict()
        for depvar, old_expr in self._solved.iteritems():
            new_expr = old_expr / sympy.diff(expr_in_old_indep, self.indep_var_symb, order)
            new_expr = new_expr.subs({self.indep_var_symb: sympy.solve(
                new_expr - new_indep_symb, self.indep_var_symb)[0]})
            new_expr = new_expr.subs({expr_in_old_indep: new_indep_symb})
            new_solved[depvar] = new_expr
        return self.__class__(self, odeqs=new_odeqs, solved=new_solved,
                              indep_var_symb=new_indep_symb)

    @property
    def is_first_order(self):
        return all([v.order == 1 for v in self._odeqs_by_dep_var.values()])

    @property
    def is_autonomous(self):
        for dep_var_func_symb, (order, expr) in \
                self._odeqs_by_dep_var.iteritems():
            unfunc_subs = {dvfs: sympy.symbol(dvfs.func.__name__) for \
                           dvfs in self.all_depvar}
            if expr.subs(unfunc_subs).diff(self.indep_var_symb) != 0:
                return False
        return True

    @property
    def is_linear(self):
        # Most easily done for first order system?
        for dep_var_func_symb, (order, expr) in \
                self._odeqs_by_dep_var.iteritems():
            for wrt in self.all_depvar:
                expr = expr.diff(wrt)

            unfunc_subs = {dvfs: sympy.symbol(dvfs.func.__name__) for \
                           dvfs in self.all_depvar}
            if expr.subs(unfunc_subs).diff(self.indep_var_symb) != 0:
                return False
        return True

    @property
    def is_homogeneous(self):
        pass

    def do_sanity_check_of_odeqs(self):
        for fnc, (order, expr) in self._odeqs_by_dep_var.iteritems():
            # Sanity check (indentation level excluded ;-)
            for check_fnc in self.all_depvar:
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
        return [self.eq(depvar) for depvar in self.all_depvar]

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
            param_symbs = param_symbs.union(get_new_symbs(
                expr, odeqs_by_dep_var.keys() + \
                [indep_var_symb]))
        new_instance = cls(odeqs_by_dep_var, indep_var_symb,
                           list(param_symbs))
        new_instance.do_sanity_check_of_odeqs()
        return new_instance


class AnyOrderODESystem(ODESystem):

    _attrs_to_cmp_for_eq = ODESystem._attrs_to_cmp_for_eq +\
                           ['_odeqs_by_dep_var']

    @staticmethod
    def mk_odeqs_by_dep_var():
        """ Convenience function for instantiating OrderedDefaultdict(ODEExpr) """
        return OrderedDefaultdict(ODEExpr)


    def __init__(self, odeqs_by_dep_var, indep_var_symb, param_symbs, solved):
        super(AnyOrderODESystem, self).__init__()
        self._odeqs_by_dep_var = odeqs_by_dep_var
        self.indep_var_symb = indep_var_symb
        self.param_symbs = param_symbs
        self._solved = solved

    def __getitem__(self, key):
        match = None
        for known_symb in self.known_symbs:
            if str(known_symb) == str(key):
                if match == None:
                    match = known_symb
                else:
                    raise KeyError('Key ambigous, there are several symbols with same str repr')
        if match == None:
            raise KeyError('Key not found: {}'.format(key))
        return match


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

    _dep_var_func_symbs = None

    _attrs_to_cmp_for_eq = ODESystem._attrs_to_cmp_for_eq +\
                           ['f', '_dep_var_func_symbs']

    is_first_order = True # no need to check that in this case

    # TODO: implement the routines for variable substitution

    def __init__(self, odesys = None):
        # Using property functions makes overwriting harder therefore
        # 'indep_var_symb', '_dep_var_func_symbs', 'param_symbs', 'f'
        # are initialized via instance methods that by default initializes
        # empty list / dict only if attribute is missing
        super(FirstOrderODESystem, self).__init__()
        if odesys != None:
            assert odesys.is_first_order
            copy_attrs = ['indep_var_symb', '_1st_ordr_red_helper_fncs',
                          'all_depvar', 'param_symbs', 'f', '_solved']
            for attr in copy_attrs:
                setattr(self, attr, getattr(odesys, attr))
        else:
            self._init_dep_var_func_symbs()
            self._init_param_symbs()
            self.init_f()

    @property
    def all_depvar(self):
        return self._dep_var_func_symbs

    @all_depvar.setter
    def all_depvar(self, value):
        self._dep_var_func_symbs = value

    @property
    def _odeqs_by_dep_var(self):
        return OrderedDefaultdict(ODEExpr, [(k, ODEExpr(1, self.f[k])) for k \
                           in self._dep_var_func_symbs])

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
                expr = expr.subs({yi: expr for yi, (expr, sol_symbs) \
                                  in self._solved.iteritems()})
                Jac_row_off_diag = [expr.diff(m) for m \
                                    in self._dep_var_func_symbs if m != yi]
                if all([c == 0 for c in Jac_row_off_diag]):
                    # Attempt solution (actually: assume success)
                    rel = sympy.Eq(yi.diff(self.indep_var_symb), expr)
                    sol = sympy.dsolve(rel, yi)
                    # Assign new symbol to inital value
                    sol_symbs = get_new_symbs(sol.rhs, self.known_symbs)
                    self._solved[yi] = sol.rhs, sol_symbs
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
        if self._dep_var_func_symbs == None:
            self._dep_var_func_symbs = []

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
        len(self._dep_var_func_symb) for the first-order derivatives
        of the self._dep_var_func_symbs (the underived dep_var_func_symb
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
        Transforms the dict mapping _dep_var_func_symbs self.f
        """
        return sympy.Matrix(
            1, len(self._dep_var_func_symbs),
            lambda q, i: self.f[self._dep_var_func_symbs[i]]
            )

    def _get_all_num_subs(self, indep_val, dep_vals, param_vals):
        indep_subs = {self.indep_var_symb: indep_val}
        dep_subs = dict(zip(self._dep_var_func_symbs, dep_vals))
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
        in self.indep_var_symb, self._dep_var_func_symbs, self._param_symbs
        for provided values of the independent, the dependent and
        parameter variables (with same order)

        Note: The signature of the function employs float point data
              (or lists thereof) in order to be compatible with
              e.g. scipy integrators, hence _get_all_numb_subs
        """
        all_subs = self._get_all_num_subs(indep_val, dep_vals, param_vals)
        return [x.subs(all_subs) for x in \
                [self.f[k] for k in self.non_analytic_depvar]]


    def dydt_jac(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating jacobian of dydt for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self._dep_var_func_symbs,
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
                self._fmat().jacobian(self.non_analytic_depvar).tolist()]


    def d2ydt2(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating d2ydt2 (f) for
        provided values which is used for substituting the the symbols
        in self.indep_var_symb, self._dep_var_func_symbs,
        self._param_symbs for provided values of the independent,
        the dependent and parameter variables (provided as dictionaries

        Note: The signature of the function employs float point
        data (or lists thereof) in order to be compatible with
        e.g. scipy integrators, hence _get_all_numb_subs
        """
        f = {k.diff(self.indep_var_symb): v for k, v in \
             self.f.iteritems()}
        dfdt_lst = [self.f[y].diff(self.indep_var_symb).subs(f) for \
                    y in self.non_analytic_depvar]
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
        self.param_tokens = self.param_tokens or []
        super(SimpleFirstOrderODESystem, self).__init__(*args, **kwargs)

    def get_param_vals_by_symb_from_by_token(self, param_vals_by_token):
        for token in param_vals_by_token:
            if not token in self.param_tokens: raise KeyError(
                'Parameter token ``{}" unknown'.format(token))
        return {self[k]: v for k, v in param_vals_by_token.iteritems()}

    def __getitem__(self, key):
        """
        If one wants to access the symbol of a _dep_var_func_symbs
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
        assert self._dep_var_func_symbs == None
        self._dep_var_func_symbs = [self[y] for y in self.dep_var_tokens]

    def _init_param_symbs(self):
        # The assert is there to signal need to subclass if using
        # SimpleFirstOrderODESystem
        assert self.param_symbs == None
        self.param_symbs = [self[k] for k in self.param_tokens]

    # End overloading
