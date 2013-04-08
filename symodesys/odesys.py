from symodesys.helpers import OrderedDefaultdict, subs_set, get_new_symbs

import sympy

from functools import reduce
from operator import add
from collections import namedtuple, OrderedDict

#This might be implemented+named better (unclear if namedtuple appropr):
ODEExpr = namedtuple('ODEExpr', 'order expr')


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


    #_attrs_to_cmp_for_eq is used for checking equality of instances
    # of ODESystem and its subclasses, the latter requires the
    # list of attrs be a common subset of attributes
    _attrs_to_cmp_for_eq = ['eqs', 'param_symbs', 'solved', 'all_depv']

    # _canonical_attrs specifies what attributes are the _defining_
    # attributes of the (sub)class. Hence it maybe changed in subclasses
    _canonical_attrs = ['odeqs', 'indepv', 'param_symbs',
                        'solved', 'frst_red_hlprs']


    def __init__(self, odesys = None, **kwargs):
        """
        Arguments:
        - `kwargs`: can be any of the keys in self._canonical_attrs
        """
        for attr in self._canonical_attrs:
            if attr in kwargs:
                setattr(self, attr, kwargs[attr])
        for attr in self._canonical_attrs:
            if not attr in kwargs:
                if odesys != None:
                    setattr(self, attr, getattr(odesys, attr))
        # Idempotent initiation of attributes
        self.solved = self.solved or {}
        self.frst_red_hlprs = self.frst_red_hlprs or []

    def __getitem__(self, key):
        if isinstance(key, sympy.Basic):
            match = None
            for known_symb in self.known_symbs:
                if str(known_symb) == str(key):
                    if match == None:
                        match = known_symb
                    else:
                        raise KeyError(
                            'Key ambigous, there are ' +\
                            'several symbols with same str repr')
            if match == None:
                raise KeyError('Key not found: {}'.format(key))
        else:
            try:
                return self[sympy.Symbol(key)]
            except KeyError:
                return self[sympy.Function(key)(self.indepv)]
        return match


    @property
    def all_depv(self):
        return self.odeqs.keys()

    @property
    def non_analytic_depv(self):
        return [x for x in self.all_depv if not x in self.solved]

    @property
    def analytic_depv(self):
        return [x for x in self.all_depv if x in self.solved]


    @property
    def solved_exprs(self):
        """ Convenience attribute """
        return [self.solved[yi][0] for yi in self.analytic_depv]

    @property
    def param_and_sol_symbs(self):
        """ Convenience attribute """
        analytic_sol_symbs = set()
        if len(self.solved) > 0:
            for expr, sol_symbs in self.solved.values():
                analytic_sol_symbs = analytic_sol_symbs.union(sol_symbs)
        analytic_sol_symbs = list(analytic_sol_symbs)
        return self.param_symbs + analytic_sol_symbs

    @property
    def known_symbs(self):
        """ Convenience attribute """
        return [self.indepv] + self.all_depv + self.param_and_sol_symbs

    def subs(self, subsd):
        for key in subsd.keys():
            if not key in self.known_symbs:
                raise KeyError(
                    'Symbol: {} not known in odesys'.format(key))

            if key in self.all_depv:
                raise KeyError('Symbol: {} is a dependent variable,' + \
                               ' use `transform_depvs`'.format(key))
            elif key == self.indepv:
                self.transform_indep(subsd[key], key)
            elif key in self.param_symbs:
                i = self.param_symbs.index(key)
                self.param_symbs[i] = self.param_symbs[i].subs(
                    {key: subsd[key]})
                for key, odeexpr in self.odeqs.iteritems():
                    self.odeqs[key] = ODEExpr(
                        odeexpr.order, odeexpr.expr.subs(
                            {key: subsd[key]}))
                for key, (expr, sol_symbs) in self.solved.iteritems():
                    self.solved[key] = expr.subs({key: subsd[key]}),\
                                       subs_set(sol_symbs,
                                                {key: subsd[key]})


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
        """ Convenience function for instantiating
        OrderedDefaultdict(ODEExpr) """
        return OrderedDefaultdict(ODEExpr)


    def attempt_analytic_sol(self, depv, hypoexpr, sol_consts):
        """
        Checks if provided analytic (similiar to
        sympy.solvers.ode.checkodesol) expr ``hypoexpr'' solves
        diff eq of ``depv'' in odesys. If it does new symbols
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

    f = None

    _attrs_to_cmp_for_eq = ODESystem._attrs_to_cmp_for_eq +\
                           ['f']

    _canonical_attrs = ['f', 'indepv', 'param_symbs',
                        'solved', 'frst_red_hlprs']


    def __init__(self, odesys = None, **kwargs):
        # Using property functions makes overwriting harder therefore
        # 'indepv', 'all_depv', 'param_symbs', 'f'
        # are initialized via instance methods that by default
        # initializes empty list / dict only if attribute is missing
        super(FirstOrderODESystem, self).__init__(odesys, **kwargs)
        self._init_param_symbs()
        if self.f == None:
            self.init_f()
        assert self.is_first_order

    @property
    def all_depv(self):
        return self.f.keys()

    @property
    def odeqs(self):
        return OrderedDefaultdict(ODEExpr, [
            (k, ODEExpr(1, self.f[k])) for k in self.all_depv])

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
                                    in self.all_depv if m != yi]
                if all([c == 0 for c in Jac_row_off_diag]):
                    # Attempt solution (actually: assume success)
                    rel = sympy.Eq(yi.diff(self.indepv), expr)
                    sol = sympy.dsolve(rel, yi)
                    # Assign new symbol to inital value
                    sol_symbs = get_new_symbs(sol.rhs, self.known_symbs)
                    self.solved[yi] = sol.rhs, sol_symbs
                    changed_last_loop = True


    def _init_param_symbs(self):
        """
        To be subclassed (or add list prop: param_symbs)

        should return list of sympy.symbols(``token_string'') instances
        The order in this list defines indices in vectors and matrices
        used by underlying numerical integration.
        (When subclassing, sympy.symarray might be useful.)
        """
        if self.param_symbs == None:
            self.param_symbs = []

    def init_f(self):
        """
        To be subclassed.

        *Is only exectuted if and only if self.f != None
        *self.init_f() must:
          set self.f to a OrderedDict with the first-order
          derivatives as values (and dependent variable sympy.Function
          instances as keys)
        """
        self.f = {}

    def _fmat(self):
        """
        Convert self.f to sympy Matrix
        """
        return sympy.Matrix(
            1, len(self.all_depv),
            lambda q, i: self.f[self.all_depv[i]]
            )

    @property
    def non_analytic_f(self):
        if len(self.solved) > 0:
            return OrderedDict(
                [(k,  self.f[k].subs(self.solved)) for k in \
                 self.non_analytic_depv])
        else:
            return self.f

    @property
    def _non_analytic_fmat(self):
        return sympy.Matrix(
            1, len(self.non_analytic_depv),
            lambda q, i: self.f[self.non_analytic_depv[i]]
            )

    @property
    def non_analytic_jac(self):
        return self._non_analytic_fmat.jacobian(self.non_analytic_depv)

    def _get_all_num_subs(self, indep_val, dep_vals, param_vals):
        indep_subs = {self.indepv: indep_val}
        dep_subs = dict(zip(self.all_depv, dep_vals))
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
        in self.indepv, self.all_depv, self._param_symbs
        for provided values of the independent, the dependent and
        parameter variables (with same order)

        Note: The signature of the function employs float point data
              (or lists thereof) in order to be compatible with
              e.g. scipy integrators, hence _get_all_numb_subs
        """
        all_subs = self._get_all_num_subs(indep_val, dep_vals,
                                          param_vals)
        return [x.subs(all_subs) for x in \
                [self.f[k] for k in self.non_analytic_depv]]


    def dydt_jac(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating jacobian of dydt for
        provided values which is used for substituting the the symbols
        in self.indepv, self.all_depv,
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
                self.non_analytic_jac.tolist()]


    def d2ydt2(self, indep_val, dep_vals, param_vals):
        """
        Convenience function for evaluating d2ydt2 (f) for
        provided values which is used for substituting the the symbols
        in self.indepv, self.all_depv,
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


    def transform_depv(self, trnsfm, inv_trnsfm):
        """
        trnsfm: dict mapping old_depv to tuple of (new_depv, expr_in_old)
        inv_subs: dict mapping old_depv to expression in new_depv
        """
        new_f = OrderedDict()
        for old_depv, old_expr in self.f.iteritems():
            if not old_depv in inv_trnsfm:
                new_f[old_depv] = self.f[old_depv].subs(inv_subs)
                continue
        for new_depv, rel in trnsfm.iteritems():
            new_rel = rel.diff(self.indepv)
            eqsd = {eq.lhs: eq.rhs for eq in self.eqs}
            new_expr = new_rel.subs(eqsd)
            new_expr = new_expr.subs(inv_trnsfm)
            new_f[new_depv] = new_expr

        # Handle solved:
        new_solved = {}
        for old_depv, (old_expr, sol_symbs) in self.solved.iteritems():
            if not old_depv in trnsfm:
                # Save old analytic expr (expressed in new
                # depv through inv_subs)
                new_solved[old_depv] = old_expr.subs(inv_subs), sol_symbs
                continue
            new_depv, trnsfm_expr = trnsfm[old_depv]
            # Save new analytic expr (expressed in new depv
            # through inv_subs)
            new_solved[new_depv] = old_expr.subs(inv_subs), sol_symbs
        # Return new instance based on carbon copy of self
        new_instance = self.__class__(self, f=new_f, solved=new_solved)
        return new_instance

    def transform_indep(self, new_indep_symb, expr_in_old_indep):
        new_f = OrderedDict()
        for depv, old_expr in self.f.iteritems():
            new_expr = old_expr / \
                       sympy.diff(expr_in_old_indep, self.indepv)
            new_expr = new_expr.subs({self.indepv: sympy.solve(
                new_expr - new_indep_symb, self.indepv)[0]})
            new_expr = new_expr.subs({expr_in_old_indep: new_indep_symb})
            new_f[depv] = new_expr

        # Handle solved:
        new_solved = {}
        for depv, (old_expr, sol_symbs) in self.solved.iteritems():
            new_expr = old_expr / sympy.diff(expr_in_old_indep,
                                             self.indepv, order)
            new_expr = new_expr.subs({self.indepv: sympy.solve(
                new_expr - new_indep_symb, self.indepv)[0]})
            new_expr = new_expr.subs({expr_in_old_indep: new_indep_symb})
            new_solved[depv] = new_expr, sol_symbs
        return self.__class__(self, f=new_f, solved=new_solved,
                              indepv=new_indep_symb)


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
    expressions = None

    def __init__(self, *args, **kwargs):
        self.param_tokens = self.param_tokens or []
        super(SimpleFirstOrderODESystem, self).__init__(*args, **kwargs)

    def get_param_vals_by_symb_from_by_token(self, param_vals_by_token):
        """
        Convenience function
        """
        for token in param_vals_by_token:
            if not token in self.param_tokens: raise KeyError(
                'Parameter token ``{}" unknown'.format(token))
        return {self[k]: v for k, v in param_vals_by_token.iteritems()}

    # Begin overloading

    def init_f(self):
        # First we need to set the keys (needed when self.expressions()
        # makes look-ups)
        self.f = OrderedDict(
            [(sympy.Function(tok)(self.indepv), None) for\
                              tok in self.dep_var_tokens])
        for tok in self.dep_var_tokens:
            self.f[self[tok]] = self.expressions[self[tok]]

    def param_vals_by_symb(self, params_by_token):
        return dict([(self[k], params_by_token[k]) for\
                     k in self.param_tokens])


    def _init_param_symbs(self):
        # The assert is there to signal need to subclass if using
        # SimpleFirstOrderODESystem
        if self.param_symbs == None:
            self.param_symbs = [sympy.Symbol(k) for k in\
                                self.param_tokens]

    # End overloading
