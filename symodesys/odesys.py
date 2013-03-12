from symodesys.helpers import OrderedDefaultdict, subs_set, get_new_symbs

import sympy

from functools import reduce
from operator import add
from collections import namedtuple

ODEExpr = namedtuple('ODEExpr', 'order expr') #<-- Might be implemented+named better


class ODESystem(object):
    """
    Central to the ODESystem if the attribute odeqs
    which should be an OrderedDefaultdict(ODEExpr). The keys in the
    dictionary should be of the form:
       sympy.Function('name')(self.indepv)
    """

    # For tracking what dep var func symbs are generated upon
    # order reduction (instantiate as list)
    # it is useful if one does not want to e.g. plot them
    frst_red_hlprs = None

    # We need holder for analytically solved parts (instantiate to dict):
    solved = None

    indepv = None # ODE implies 1 indep. variable,
                          #  set to sympy.symbol(...)
    param_symbs = None
    f = None


    #_attrs_to_cmp_for_eq is used for checking equality of instances of ODESystem
    # and its subclasses, the latter requires the list of attrs be a common subset
    # of attributes
    _attrs_to_cmp_for_eq = ['eqs', 'param_symbs', 'solved', 'all_depv']

    # _canonical_attrs specifies what attributes are the _defining_ attributes of the
    # (sub)class. Hence it maybe changed in subclasses
    _canonical_attrs = ['odeqs', 'indepv', 'param_symbs', 'solved', 'frst_red_hlprs']


    def __init__(self, odesys = None, **kwargs):
        """
        Arguments:
        - `kwargs`: can be any of the keys in self._canonical_attrs
        """
        for attr in self._canonical_attrs:
            if attr in kwargs:
                setattr(self, attr, kwargs[attr])
            else:
                if odesys != None:
                    setattr(self, attr, getattr(odesys, attr))
        # Idempotent initiation of attributes
        self.solved = self.solved or {}
        self.frst_red_hlprs = self.frst_red_hlprs or []


    @property
    def non_analytic_depv(self):
        return [x for x in self.odeqs if not x in self.solved]

    @property
    def analytic_depv(self):
        return [x for x in self.odeqs if x in self.solved]

    @property
    def all_depv(self):
        return self.odeqs.keys()

    @property
    def solved_exprs(self):
        """ Convenience attribute """
        return [self.solved[yi][0] for yi in self.analytic_depv]

    @property
    def known_symbs(self):
        """ Convenience attribute """
        analytic_sol_symbs = set()
        if len(self.solved) > 0:
            for expr, sol_symbs in self.solved.values():
                analytic_sol_symbs = analytic_sol_symbs.union(sol_symbs)
        analytic_sol_symbs = list(analytic_sol_symbs)
        return [self.indepv] + self.all_depv +\
               self.param_symbs + analytic_sol_symbs

    def subs(self, subsd):
        # TODO: fix special cases self.indepv and self.depv
        for key in subsd.keys():
            if not key in self.known_symbs:
                raise KeyError('Symbol: {} not known in odesys'.format(key))

            if key in self.all_depv:
                raise KeyError('Symbol: {} is a dependent variable,' + \
                               ' use `transform_depvs`'.format(key))
            elif key == self.indepv:
                self.transform_indep(subsd[key], key)
            elif key in self.param_symbs:
                i = self.param_symbs.index(key)
                self.param_symbs[i] = self.param_symbs[i].subs({key: subsd[key]})
                for key, odeexpr in self.odeqs.iteritems():
                    self.odeqs[key] = ODEExpr(
                        odeexpr.order, odeexpr.expr.subs({key: subsd[key]}))
                for key, (expr, sol_symbs) in self.solved.iteritems():
                    self.solved[key] = expr.subs({key: subsd[key]}),\
                                        subs_set(sol_symbs, {key: subsd[key]})


    def transform_depv(self, trnsf, inv_subs):
        """
        trnsf: dict mapping old_depv to tuple of (new_depv, expr_in_old)
        inv_subs: dict mapping old_depv to expression in new_depv
        """
        new_odeqs = OrderedDefaultdict(ODEExpr)
        for old_depv, (order, old_expr) in self.odeqs.iteritems():
            if not old_depv in trnsf:
                new_odeqs[old_depv] = self.odeqs[old_depv].subs(inv_subs)
                continue
            # Unpack forward transform
            new_depv, trnsf_expr = trnsf[old_depv]
            eqsd = {eq.lhs: eq.rhs for eq in self.eqs}
            # Get new diff eq
            diff_expr = trnsf_expr.diff(self.indepv, order).subs(eqsd)
            # Save new diff eq (expressed in new depv through inv_subs)
            new_odeqs[new_depv] = ODEExpr(order, diff_expr.subs(inv_subs))

        # Handle solved:
        new_solved = {}
        for old_depv, (old_expr, sol_symbs) in self.solved.iteritems():
            if not old_depv in trnsf:
                # Save old analytic expr (expressed in new depv through inv_subs)
                new_solved[old_depv] = old_expr.subs(inv_subs), sol_symbs
                continue
            new_depv, trnsf_expr = trnsf[old_depv]
            # Save new analytic expr (expressed in new depv through inv_subs)
            new_solved[new_depv] = old_expr.subs(inv_subs), sol_symbs
        # Return new instance based on carbon copy of self
        return self.__class__(self, odeqs=new_odeqs, solved=new_solved)

    def transform_indep(self, new_indep_symb, expr_in_old_indep):
        new_odeqs = OrderedDefaultdict(ODEExpr)
        for depv, (order, old_expr) in self.odeqs.iteritems():
            new_expr = old_expr / sympy.diff(expr_in_old_indep, self.indepv, order)
            new_expr = new_expr.subs({self.indepv: sympy.solve(
                new_expr - new_indep_symb, self.indepv)[0]})
            new_expr = new_expr.subs({expr_in_old_indep: new_indep_symb})
            new_odeqs[depv] = ODEExpr(order, new_expr)

        # Handle solved:
        new_solved = {}
        for depv, (old_expr, sol_symbs) in self.solved.iteritems():
            new_expr = old_expr / sympy.diff(expr_in_old_indep, self.indepv, order)
            new_expr = new_expr.subs({self.indepv: sympy.solve(
                new_expr - new_indep_symb, self.indepv)[0]})
            new_expr = new_expr.subs({expr_in_old_indep: new_indep_symb})
            new_solved[depv] = new_expr, sol_symbs
        return self.__class__(self, odeqs=new_odeqs, solved=new_solved,
                              indepv=new_indep_symb)

    @property
    def is_first_order(self):
        return all([v.order == 1 for v in self.odeqs.values()])

    @property
    def is_autonomous(self):
        for depv, (order, expr) in \
                self.odeqs.iteritems():
            unfunc_subs = {dvfs: sympy.symbol(dvfs.func.__name__) for \
                           dvfs in self.all_depv}
            if expr.subs(unfunc_subs).diff(self.indepv) != 0:
                return False
        return True

    @property
    def is_linear(self):
        # Most easily done for first order system?
        for depv, (order, expr) in \
                self.odeqs.iteritems():
            for wrt in self.all_depv:
                expr = expr.diff(wrt)

            unfunc_subs = {dvfs: sympy.symbol(dvfs.func.__name__) for \
                           dvfs in self.all_depv}
            if expr.subs(unfunc_subs).diff(self.indepv) != 0:
                return False
        return True

    @property
    def is_homogeneous(self):
        pass

    def do_sanity_check_of_odeqs(self):
        for fnc, (order, expr) in self.odeqs.iteritems():
            # Sanity check (indentation level excluded ;-)
            for check_fnc in self.all_depv:
                # Make sure expr don't contain derivatives larger
                # than in odeqs order:
                if expr.has(check_fnc):
                    for arg in expr.args:
                        if arg.has(check_fnc):
                            if arg.is_Derivative:
                                fnc, wrt = args[0], args[1:]
                                assert not len(wrt) > \
                                       odeqs[check_fnc][0]

    def __eq__(self, other):
        for attr in self._attrs_to_cmp_for_eq:
            if getattr(self, attr) != getattr(other, attr): return False
        return True

    @property
    def eqs(self):
        """
        Returns a list of Sympy Eq instances describing the ODE system
        """
        return [self.eq(depv) for depv in self.all_depv]

    def eq(self, depv):
        """
        Returns a sympy.Eq for diff eq of ``depv''
        """
        order, expr = self.odeqs[depv]
        return sympy.Eq(depv.diff(self.indepv, order), expr)

    @classmethod
    def from_list_of_eqs(cls, lst):
        """
        Determine independnt variable, dependent variables
        and parameter variables and perform Sanity checks
        """
        fncs = []

        # independet variable as key, (order, expr) as value
        odeqs = OrderedDefaultdict(ODEExpr)
        indepv = None
        for eq in lst:
            assert eq.lhs.is_Derivative
            fnc, wrt = eq.lhs.args[0], eq.lhs.args[1:]
            assert fnc not in odeqs
            if indepv == None:
                assert all([wrt[0] == x for x in wrt]) # No PDEs!
                indepv = wrt[0]
            else:
                assert all([indepv == x for x in wrt])
            odeqs[fnc] = ODEExpr(len(wrt), eq.rhs)

        param_symbs = set()
        for order, expr in odeqs.values():
            param_symbs = param_symbs.union(get_new_symbs(
                expr, odeqs.keys() + \
                [indepv]))
        new_instance = cls(odeqs, indepv,
                           list(param_symbs))
        new_instance.do_sanity_check_of_odeqs()
        return new_instance


class AnyOrderODESystem(ODESystem):

    _attrs_to_cmp_for_eq = ODESystem._attrs_to_cmp_for_eq +\
                           ['odeqs']

    @staticmethod
    def mkodeqs():
        """ Convenience function for instantiating OrderedDefaultdict(ODEExpr) """
        return OrderedDefaultdict(ODEExpr)

    def __getitem__(self, key):
        if isinstance(key, sympy.Basic):
            match = None
            for known_symb in self.known_symbs:
                if str(known_symb) == str(key):
                    if match == None:
                        match = known_symb
                    else:
                        raise KeyError('Key ambigous, there are several symbols with same str repr')
            if match == None:
                raise KeyError('Key not found: {}'.format(key))
        else:
            return self[sympy.Symbol(key)]
        return match


    def attempt_analytic_sol(self, depv, hypoexpr, sol_consts):
        """
        Checks if provided analytic (similiar to sympy.solvers.ode.checkodesol)
        expr ``hypoexpr'' solves diff eq of ``depv'' in odesys. If it does
        new symbols
        """
        eq = self.eq(depv).subs({depv: hypoexpr})
        if bool(eq.doit()):
            hypoparams = get_new_symbs(hypoexpr, self.known_symbs)
            self.solved[depv] = hypoexpr, sol_consts
            return True
        else:
            return False


    @property
    def f(self):
        assert self.is_first_order
        return OrderedDefaultdict(ODEExpr, [\
            (k, v.expr) for k, v in self.odeqs.items()])

    def get_helper_fnc(self, fnc, order):
        """
        Returns a list of of length order - 1
        for use in reformulation of higher order
        ODE eq in first order ODE (sub) sys.
        """
        helpers = {}
        for i in range(1, order):
            candidate = sympy.Function(str(fnc.func.__name__) + '_h' + \
                                       str(i))(self.indepv)
            while candidate in self.odeqs:
                candidate = candidate + '_h'
            helpers[i] = candidate
        return helpers

    def reduce_to_sys_of_first_order(self):
        """ Returns a new instance with reduced order (1st) """
        new_odeqs = self.mkodeqs()
        frst_red_hlprs = self.frst_red_hlprs[:]
        # TODO, revise frst_red_hlprs format (cur.  len 3 tuple)
        # and update analytic_harmonic_oscillator.py and IVP.plot()
        for fnc, (order, expr) in self.odeqs.iteritems():
            if order == 1:
                new_odeqs[fnc] = ODEExpr(order, expr)
                continue
            hlpr = self.get_helper_fnc(fnc, order)
            new_odeqs[fnc] = ODEExpr(1, hlpr[1])
            for o in range(1, order - 1):
               new_odeqs[hlpr[o]] = ODEExpr(1, hlpr[o + 1])
            subsd = {sympy.Derivative(fnc, self.indepv, i): \
                     hlpr[i] for i in range(1, order)}
            new_odeqs[hlpr[order - 1]] = ODEExpr(1, expr.subs(subsd))
            frst_red_hlprs.extend([
                (fnc, o, expr) for o, expr in hlpr.iteritems()])
        new_instance = self.__class__(new_odeqs, self.indepv,
                                      self.param_symbs[:])
        new_instance.frst_red_hlprs = frst_red_hlprs
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

    depv = None

    _attrs_to_cmp_for_eq = ODESystem._attrs_to_cmp_for_eq +\
                           ['f', 'depv']

    _canonical_attrs = ['f', 'indepv', 'param_symbs', 'solved', 'frst_red_hlprs']

    # TODO: implement the routines for variable substitution
    # TODO: Require `f` to be OrderedDict

    def __init__(self, odesys = None, **kwargs):
        # Using property functions makes overwriting harder therefore
        # 'indepv', 'depv', 'param_symbs', 'f'
        # are initialized via instance methods that by default initializes
        # empty list / dict only if attribute is missing
        super(FirstOrderODESystem, self).__init__(odesys, **kwargs)
        self._init_depv()
        self._init_param_symbs()
        self.init_f()
        assert self.is_first_order

    @property
    def all_depv(self):
        return self.depv

    # @all_depv.setter
    # def all_depv(self, value):
    #     self.depv = value

    @property
    def odeqs(self):
        return OrderedDefaultdict(ODEExpr, [(k, ODEExpr(1, self.f[k])) for k \
                           in self.depv])

    def recursive_analytic_auto_sol(self):
        """
        Solves equations one by one
        """
        changed_last_loop = True
        # new_y = []
        while changed_last_loop:
            changed_last_loop = False
            for yi, expr in self.f.iteritems():
                if yi in self.solved: continue
                # Check to see if it only depends on itself
                expr = expr.subs({yi: expr for yi, (expr, sol_symbs) \
                                  in self.solved.iteritems()})
                Jac_row_off_diag = [expr.diff(m) for m \
                                    in self.depv if m != yi]
                if all([c == 0 for c in Jac_row_off_diag]):
                    # Attempt solution (actually: assume success)
                    rel = sympy.Eq(yi.diff(self.indepv), expr)
                    sol = sympy.dsolve(rel, yi)
                    # Assign new symbol to inital value
                    sol_symbs = get_new_symbs(sol.rhs, self.known_symbs)
                    self.solved[yi] = sol.rhs, sol_symbs
                    changed_last_loop = True
                    # new_y.append(yi)
        # return new_y

    def _init_depv(self):
        """
        To be subclassed (or add list prop: dep_var_func_symbs)

        should return list of
          sympy.Function(``token_string'')(self.indep_var) instances
        The order in this list defines indices in vectors and matrices
        used by underlying numerical integration.
        """
        if self.depv == None:
            self.depv = []

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

        self.f should return a OrderedDict with the first-order
        derivatives as values (and dependent sympy name as)
        of the self.depv (the underived dep_var_func_symb
        acts as key) expressed solely in numerical constants, sympy
        function expressions, indepv, dep_var_func_symbs and
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
        Transforms the dict mapping depv self.f
        """
        return sympy.Matrix(
            1, len(self.depv),
            lambda q, i: self.f[self.depv[i]]
            )

    def _get_all_num_subs(self, indep_val, dep_vals, param_vals):
        indep_subs = {self.indepv: indep_val}
        dep_subs = dict(zip(self.depv, dep_vals))
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
        in self.indepv, self.depv, self._param_symbs
        for provided values of the independent, the dependent and
        parameter variables (with same order)

        Note: The signature of the function employs float point data
              (or lists thereof) in order to be compatible with
              e.g. scipy integrators, hence _get_all_numb_subs
        """
        all_subs = self._get_all_num_subs(indep_val, dep_vals, param_vals)
        return [x.subs(all_subs) for x in \
                [self.f[k] for k in self.non_analytic_depv]]


    def dydt_jac(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating jacobian of dydt for
        provided values which is used for substituting the the symbols
        in self.indepv, self.depv,
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
                self._fmat().jacobian(self.non_analytic_depv).tolist()]


    def d2ydt2(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating d2ydt2 (f) for
        provided values which is used for substituting the the symbols
        in self.indepv, self.depv,
        self._param_symbs for provided values of the independent,
        the dependent and parameter variables (provided as dictionaries

        Note: The signature of the function employs float point
        data (or lists thereof) in order to be compatible with
        e.g. scipy integrators, hence _get_all_numb_subs
        """
        f = {k.diff(self.indepv): v for k, v in \
             self.f.iteritems()}
        dfdt_lst = [self.f[y].diff(self.indepv).subs(f) for \
                    y in self.non_analytic_depv]
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

    indepv = sympy.symbols('t') # ODE implies 1 indep. variable

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
        If one wants to access the symbol of a depv
        or a param_symbs and do not want to hardcode the order in
        the code for item access, it can be retrieved using this function
        """
        if key in self.dep_var_tokens:
            assert key not in self.param_tokens
            return sympy.Function(key)(self.indepv)
        elif key in self.param_tokens:
            return sympy.symbols(key)
        else:
            raise KeyError('Unknown token')


    # Begin overloading

    def param_vals_by_symb(self, params_by_token):
        return dict([(self[k], params_by_token[k]) for\
                     k in self.param_tokens])

    def _init_depv(self):
        # The assert is there to signal need to subclass if using
        # SimpleFirstOrderODESystem
        if self.depv == None:
            self.depv = [self[y] for y in self.dep_var_tokens]

    def _init_param_symbs(self):
        # The assert is there to signal need to subclass if using
        # SimpleFirstOrderODESystem
        if self.param_symbs == None:
            self.param_symbs = [self[k] for k in self.param_tokens]

    # End overloading
