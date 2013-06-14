# stdlib imports
from functools import reduce
from operator import add
from collections import namedtuple, OrderedDict
import new
import logging

# other imports
import sympy

# project imports
from symodesys.helpers import subs_set, get_new_symbs, deprecated

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def log_call_debug(func):
    from functools import wraps
    @wraps(func)
    def wrapper(*args, **kwargs):
        logging.debug(func.__name__ + 'called.')
        return func(*args, **kwargs)
    return wrapper


#This might be implemented+named better (unclear if namedtuple appropriate/adds value):
ODEExpr = namedtuple('ODEExpr', 'order expr')


class _ODESystemBase(object):
    """
    Central to the _ODESystemBase if the attribute _odeqs
    which should be an OrderedDict(). The keys in the
    dictionary should be of the form:
       sympy.Function('name')(self.indepv)

    User's should use or subclass either FirstOrderODESystem
    or AnyOrderODESystem. _ODESystemBase defines common attributes
    and methods of those derived classes.
    """

    # By default Symbols are taken to represent real valued
    # variables, override this by changing `real` to False:
    real = True

    # For tracking what dep var func symbs are generated upon
    # order reduction (instantiate as list)
    # it is useful if one does not want to e.g. plot them
    frst_red_hlprs = None

    # We need holder for analytically solved parts (instantiate to dict):
    _solved = None

    indepv = None # ODE implies 1 indep. variable,
                          #  set to sympy.Symbol(...)
    param_symbs = None


    #_attrs_to_cmp_for_eq is used for checking equality of instances
    # of _ODESystemBase and its subclasses, the latter requires the
    # list of attrs be a common subset of attributes
    _attrs_to_cmp_for_eq = ['_odeqs', 'param_symbs', '_solved', 'all_depv']

    # _canonical_attrs specifies what attributes are the _defining_
    # attributes of the (sub)class. Hence it maybe changed in subclasses
    _canonical_attrs = ['_odeqs', 'indepv', 'param_symbs',
                        '_solved', 'frst_red_hlprs'] # Minimum variables needed

    _redundant_attrs = [] # May be passed as kwargs to __init__ but must not be

    def __init__(self, odesys = None, **kwargs):
        """
        Arguments:
        - `kwargs`: can be any of the keys in self._canonical_attrs
        """
        for attr in self._canonical_attrs + self._redundant_attrs:
            if attr in kwargs:
                setattr(self, attr, kwargs[attr])
        for attr in self._canonical_attrs:
            if not attr in kwargs:
                if odesys != None:
                    setattr(self, attr, getattr(odesys, attr))
        # Idempotent initiation of attributes
        self._solved = self._solved or {}
        self._solved_undefined = []
        self.frst_red_hlprs = self.frst_red_hlprs or []

    def __getitem__(self, key):
        if isinstance(key, sympy.Basic):
            match = None
            for known_symb in self.known_symbs:
                if str(known_symb) == str(key):
                    if match == None:
                        match = known_symb
                    else:
                        # This place should never be reached
                        raise KeyError(
                            'Key ambigous, there are ' +\
                            'several symbols with same str repr')
            if match == None:
                raise KeyError('Key not found: {}'.format(key))
        else:
            try:
                return self[sympy.Symbol(key, real=self.real)]
            except KeyError:
                return self[sympy.Function(key)(self.indepv)]
        return match

    def mk_func(self, key):
        """
        Returns an sympy.Function instance with name key and correct
        dependent variable
        """
        # metaclasses and special behaviour of subclassed
        # classes of sympy.Function makes this tricky.. see tests
        try:
            ori_eval_is_real = sympy.Function._eval_is_real
        except AttributeError:
            ori_eval_is_real = None
        setattr(sympy.Function, '_eval_is_real', lambda self_: self.real)
        instance = sympy.Function(key)(self.indepv)
        if ori_eval_is_real:
            setattr(sympy.Function, '_eval_is_real', ori_eval_is_real)
        else:
            delattr(sympy.Function, '_eval_is_real')
        return instance


    def ensure_dictkeys_as_symbs(self, val_by_token):
        """
        Convenience function for converting dicts with keys of form
        'y1', 'y2' into sympy.Symbol('y1') through __getitem__ (which
        also checks existance of y1, y2... etc.)
        """
        return {self[k]: v for k, v in val_by_token.items()}


    @property
    def all_depv(self):
        """ Returns all dependent variables of the system """
        return self._odeqs.keys()


    @property
    def non_analytic_depv(self):
        """
        Returns all dependent variables of the system
        that have not been solved analytically
        """
        return [x for x in self.all_depv if not x in self._solved]


    @property
    def analytic_depv(self):
        """
        Returns all dependent variables of the system
        that have been solved analytically
        """
        return [x for x in self.all_depv if x in self._solved]


    @property
    def analytic_relations(self):
        """
        Convenience attribute, returns all the expressions of the analytic
        solutions corresponding to the dependent variables which have been solved
        analytically.
        """
        return [self._solved[yi][0] for yi in self.analytic_depv]


    @property
    def analytic_sol_symbs(self):
        """
        Returns a list of symbols introduced in when process of analytically
        solving the expressions of the dependent variables.
        """
        symbs = set()
        if len(self._solved) > 0:
            for expr, sol_symbs in self._solved.values():
                symbs = symbs.union(sol_symbs)
        return list(symbs)


    @property
    def param_and_sol_symbs(self):
        """
        Convenience attribute, this is useful when interfacing with external
        software solving ODE's since the notion of parameters is changed.
        (New parameters might have been introduced when solving some of the
        expressions of for the dependent variables analytically).
        """
        return self.param_symbs + self.analytic_sol_symbs


    @property
    def known_symbs(self):
        """
        Convenience attribute, returns a list of all Symbol/Function isntances
        in use in the system.
        """
        return [self.indepv] + self.all_depv + self.param_and_sol_symbs

    @property
    def forbidden_symbs(self):
        """
        Extends self.known_symbs to symbols with names coinciding with e.g. dependent variables.
        """
        return self.known_symbs + [sympy.Symbol(dvfs.func.__name__, real=self.real) for \
                                   dvfs in self.all_depv]


    def subs(self, subsd):
        """
        Performs variable substituions in the system according to the provided
        dictionary `subsd`
        """
        for key, value in subsd.items():
            if not key in self.known_symbs:
                raise KeyError(
                    'Symbol: {} not known in odesys'.format(key))

            if key in self.all_depv:
                raise RuntimeError(('Symbol: {} is a dependent variable,' + \
                                   ' use `transform_depvs`').format(key))
            elif key == self.indepv:
                self.transform_indep(value, key)
            elif key in self.param_and_sol_symbs:
                if key in self.param_symbs:
                    self.param_symbs[self.param_symbs.index(key)] = value
                else:
                    # Sanity check (param_and_sol_symbs)
                    assert key in self.analytic_sol_symbs
                    new_solved = {}
                    for depv, (expr, sol_symb) in self._solved.items():
                        if expr.has(key):
                            new_solved[depv] = (expr.subs({key: value}),
                                               subs_set(sol_symb, {key: value}))
                        else:
                            new_solved[depv] = (expr, sol_symb)
                    self._solved = new_solved
                for key, odeexpr in self._odeqs.iteritems():
                    self._odeqs[key] = ODEExpr(
                        odeexpr.order, odeexpr.expr.subs(
                            {key: value}))
                for key, (expr, sol_symbs) in self._solved.iteritems():
                    self._solved[key] = expr.subs({key: value}),\
                                       subs_set(sol_symbs,
                                                {key: value})
            else:
                raise KeyError("{} not known in ODE system".format(key))


    @property
    def is_first_order(self):
        """
        Returns true if the highest order of the ODE's in the system is 1.
        """
        return all([v.order == 1 for v in self._odeqs.values()])


    def get_highest_order(self):
        """
        Returns the (highest) order of the ODE System.
        If it has been completely solved analytically 0 is returned
        """
        if len(self._odeqs) > 0:
            return max((order for order, expr in self._odeqs.values()))
        else:
            return 0 # Completely solved analytically


    @property
    def is_autonomous(self):
        """
        Returns true if the system is autonomous (the independent variable is abscent
        in all expressions for the derivatives of the dependent variables)
        """
        for depv, (order, expr) in \
                self._odeqs.iteritems():
            if self.unfunc_depv(expr).diff(self.indepv) != 0:
                return False
        return True


    def unfunc_depv(self, expr):
        """
        Convenience method, transforms the Function instances of the dependent variables
        possibly present in provided argument `expr` into Symbol instances.
        """
        unfunc_subs = {dvfs: sympy.Symbol(dvfs.func.__name__, real=self.real) for \
                       dvfs in self.all_depv}
        return expr.subs(unfunc_subs)


    @property
    def is_homogeneous(self):
        """ TODO implement this """
        pass


    @log_call_debug
    def _do_sanity_check_of_odeqs(self):
        """
        Asserts that no expr contain dependent variable derivatives of higher order
        than the one currently explicitly defined for that dependent variable.
        """
        for fnc, (order, expr) in self._odeqs.iteritems():
            # Sanity check (indentation level excluded ;-)
            for check_fnc in self.all_depv:
                # Make sure expr don't contain derivatives larger
                # than in _odeqs order:
                if expr.has(check_fnc):
                    for arg in expr.args:
                        if arg.has(check_fnc):
                            if arg.is_Derivative:
                                fnc, wrt = args[0], args[1:]
                                assert not len(wrt) > \
                                       _odeqs[check_fnc][0]


    def __eq__(self, other):
        """
        Some subclasses have different canonical attributes that
        needs to be compared. These are collected in
        self._attrs_to_cmp_for_eq
        """
        for attr in self._attrs_to_cmp_for_eq:
            if getattr(self, attr) != getattr(other, attr): return False
        return True


    @property
    def eqs(self):
        """
        Returns a list of Sympy Eq instances describing the ODE system.
        It may be useful for exporting to LaTeX etc.
        """
        return [self.eq(depv) for depv in self.all_depv]


    def eq(self, depv):
        """
        Returns a sympy.Eq for diff eq of ``depv''
        """
        order, expr = self._odeqs[depv]
        return sympy.Eq(depv.diff(self.indepv, order), expr)


    @classmethod
    def from_list_of_eqs(cls, lst):
        """
        Classmethod producing an instance from a list of
        sympy Eq instances. The independent variable,
        dependent variables and parameter variables are
        determined from the Eq instances and finally
        a sanity check is performed.
        """
        fncs = []

        # independet variable as key, (order, expr) as value
        _odeqs = OrderedDict()
        indepv = None
        for eq in lst:
            assert eq.lhs.is_Derivative
            fnc, wrt = eq.lhs.args[0], eq.lhs.args[1:]
            assert fnc not in _odeqs # we cannot have multiple definitions
            if indepv == None:
                assert all([wrt[0] == x for x in wrt]) # No PDEs!
                indepv = wrt[0]
            else:
                assert all([indepv == x for x in wrt]) # System of same indepv
            _odeqs[fnc] = ODEExpr(len(wrt), eq.rhs)

        param_symbs = set()
        for order, expr in _odeqs.values():
            param_symbs = param_symbs.union(get_new_symbs(
                expr, _odeqs.keys() + \
                [indepv]))
        new_instance = cls(_odeqs=_odeqs, indepv=indepv,
                           param_symbs=list(param_symbs))
        new_instance._do_sanity_check_of_odeqs()
        return new_instance


    @log_call_debug
    def attempt_analytic_sol(self, depv, hypoexpr, sol_consts):
        """
        Checks if provided analytic (similiar to
        sympy.solvers.ode.checkodesol) expr ``hypoexpr'' solves
        diff eq of ``depv'' in odesys. If it does new symbols are
        generated and saved toghether with hypoexpr in self._solved
        """
        eq = self.eq(depv).subs({depv: hypoexpr})
        if bool(eq.doit()):
            hypoparams = get_new_symbs(hypoexpr, self.known_symbs)
            self._solved[depv] = hypoexpr, sol_consts
            return True
        else:
            return False


class AnyOrderODESystem(_ODESystemBase):
    """
    AnyOrderODESystem is a useful class that is usually
    used as an original representation which is processed
    before treated analytically or numerically (or a combiation
    thereof).
    """

    _redundant_attrs = ['f']

    @property
    def f(self):
        """
        The `f` attribute only makes sense if the System is first order
        """
        assert self.is_first_order
        return OrderedDict([\
            ODEExpr(k, v.expr) for k, v in self._odeqs.items()])

    @f.setter
    def f(self, value):
        """
        The `f` attribute assumes expressions are for
        the first order derivatives.
        """
        self._odeqs = OrderedDict()
        for depv, expr in value.items():
            self._odeqs[depv] = ODEExpr(1, expr)


    def _get_helper_fnc(self, fnc, order):
        """
        Returns a list of of length order - 1
        for use in reformulation of higher order
        ODE eq in first order ODE (sub) sys.
        """
        helpers = {}
        for i in range(1, order):
            candidate = self.mk_func(fnc.func.__name__ + '_h' + str(i))
            while candidate in self._odeqs:
                candidate = candidate + '_h'
            helpers[i] = candidate
        return helpers


    def reduce_to_sys_of_first_order(self, y0=None, default_red_init_val=None):
        """
        Returns a new instance with reduced order (1st).
        If y0 and default_red_init_val are provided, y0 will be updated
        with key,value entries (helper_function, default_red_init_val)
        """
        new_odeqs = OrderedDict()
        frst_red_hlprs = self.frst_red_hlprs[:]
        # TODO, revise frst_red_hlprs format (cur.  len 3 tuple)
        # and update analytic_harmonic_oscillator.py and IVP.plot()
        for fnc, (order, expr) in self._odeqs.iteritems():
            if order == 1:
                new_odeqs[fnc] = ODEExpr(order, expr)
                continue
            hlpr = self._get_helper_fnc(fnc, order)
            new_odeqs[fnc] = ODEExpr(1, hlpr[1])
            for o in range(1, order - 1):
               new_odeqs[hlpr[o]] = ODEExpr(1, hlpr[o + 1])
            subsd = {sympy.Derivative(fnc, self.indepv, i): \
                     hlpr[i] for i in range(1, order)}
            new_odeqs[hlpr[order - 1]] = ODEExpr(1, expr.subs(subsd))
            frst_red_hlprs.extend([
                (fnc, o, expr) for o, expr in hlpr.iteritems()])
        new_instance = FirstOrderODESystem(
            f=OrderedDict([(k, v.expr)for k,v in new_odeqs.items()]),
            indepv=self.indepv, param_symbs=self.param_symbs[:])
        new_instance.frst_red_hlprs = frst_red_hlprs
        if y0:
            for source, order, helper in frst_red_hlprs:
                if not helper in y0:
                    y0[helper] = default_red_init_val
        return new_instance


class FirstOrderODESystem(_ODESystemBase):
    """
    This class provides an efficient interface to the manipulation
    and export of first order ODE Systems, it is a special case in the
    sense that it only contains first order derivatives (which is a
    requirement for the numerical treatment of ODE Systems).

    The attribute `f` is given a special role. It is a simpliciation of _odeqs
    attribute which now has a redundant order part of its (order, expr) tuples.
    Furthermore
    Hence, `f` is a dictionary mapping the Symbol of the dependent variable to
    the expression of the first derivative of same dependent variable with respect
    to the independnt variable.

    When the ODE systems equations is generated from user
    data this class is to be subclassed to provide similar routines
    of the SimpleFirstOrderODESystem which is to be used when explicitly
    defining an ODE System.
    """

    f = None

    _canonical_attrs = ['f', 'indepv', 'param_symbs',
                        '_solved', 'frst_red_hlprs']

    _redundant_attrs = ['_odeqs']

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
    def is_linear(self):
        """
        is True if the dependent variables present in the expressions
        are at most polynomial of degree one.
        """
        for depv, (order, expr) in \
                self._odeqs.iteritems():
            for wrt in self.all_depv:
                expr = expr.diff(wrt)

            if self.unfunc_depv(expr).diff(self.indepv) != 0:
                return False
        return True


    @property
    def all_depv(self):
        """ Returns all dependent variables of the system """
        return self.f.keys()


    @property
    def _odeqs(self):
        """
        In FirstOrderODESystem f is the canonical variable and
        _odeqs can easily be generated from f.
        """
        return OrderedDict([
            ODEExpr(k, ODEExpr(1, self.f[k])) for k in self.all_depv])

    @_odeqs.setter
    def _odeqs(self, value):
        """
        In FirstOrderODESystem f is the canonical variable and
        _odeqs can easily be generated from f.
        """
        self.f = OrderedDict()
        for depv, odeexpr in value.items():
            assert odeexpr.order == 1 # must be first order
            self.f[depv] = odeexpr.expr


    @log_call_debug
    def recursive_analytic_auto_sol(self):
        """
        Solves equations one by one
        (hence it can only find solutions for
         independent equations)

        TODO: sometimes the different solutions are valid for different
        conditions of choosen parameters. Best way avoid this problem
        would be to run substitute all parameter symbols for chosen numeric
        values (see self.subs) with the obvious caveat that the ODE System
        is significantly less flexible from that point on.
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
                                    in self.all_depv if m != yi]
                if all([c == 0 for c in Jac_row_off_diag]):
                    # Attempt solution (actually: assume success)
                    rel = sympy.Eq(yi.diff(self.indepv), expr)
                    sol = sympy.dsolve(rel, yi)
                    # If sol contains a Piecewise definition,
                    # accept the default solution and store
                    # the others as undefined cases.
                    if sol.has(sympy.Piecewise):
                        args = sol.rhs.args
                        for arg in args:
                            if isinstance(arg, sympy.Piecewise):
                                # found it
                                for expr, cond in arg:
                                    if isinstance(cond, sympy.Dummy):
                                        # Someone with deeper insight
                                        # might want to improve this..
                                        assert cond.name == 'True'
                                        sol_expr = sol.rhs.fromiter(
                                            (x if x != arg else expr for x in args))
                                    else:
                                        self._solved_undefined.append(cond)
                        else:
                            raise RuntimeError("Piecewise extraction failed, improve code!")
                    else:
                        sol_expr = sol.rhs
                    # Assign new symbol to inital value
                    sol_symbs = get_new_symbs(sol_expr, self.known_symbs)
                    self._solved[yi] = sol_expr, sol_symbs
                    changed_last_loop = True


    def _init_param_symbs(self):
        """
        To be subclassed (or add attribute (list): param_symbs)

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
        self.f = OrderedDict()


    def _fmat(self):
        """
        Convert self.f to sympy Matrix, used to generate Jacobian.
        """
        return sympy.Matrix(
            1, len(self.all_depv),
            lambda q, i: self.f[self.all_depv[i]]
            )


    @property
    def non_analytic_f(self):
        if len(self._solved) > 0:
            return OrderedDict(
                [(k,  self.f[k].subs(self._solved)) for k in \
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


    def is_stiff(self, indep_val, dep_vals, param_vals, criteria=100):
        pass
        # eigenvalues = scipy.eigenvalues(dfdy)
        # mineigen = np.min(eigenvalues)
        # tspan = tend-t0
        # if mineigen < 0:
        #     if -mineigen*tspan > criteria:
        #         return True
        # return False


    def transform_depv(self, trnsfm, inv_trnsfm):
        """
        trnsfm: dict mapping new_depv to expr_in_old
        inv_trnsfm: dict mapping old_depv to expression in new_depv
        """
        new_f = OrderedDict()
        # User must give transformation for same number of variables
        assert len(trnsfm) == len(inv_trnsfm) == len(self.all_depv)
        assert all([old_depv in self.all_depv for old_depv in inv_trnsfm.keys()])

        # We do not accept that a dependent variable is lost in the transformation:
        for old_depv in self.all_depv:
            used = False
            for rel in trnsfm.values():
                if rel.has(old_depv):
                    used = True
                    break
            assert used

        eqsd = {eq.lhs: eq.rhs for eq in self.eqs} # d/dt(old): expr_in_old
        new_solved = {}
        for new_depv, rel_in_old in trnsfm.iteritems():
            # First use our expressions for solved variables
            rel_in_old = rel_in_old.subs(dict(zip(
                self.analytic_depv, self.analytic_relations)))
            expr_in_old = rel_in_old.diff(self.indepv) # what does the derivatives look like?
            analytic = True
            for depv in self.non_analytic_depv:
                if expr_in_old.has(depv.diff(self.indepv)):
                    analytic = False
                    break
            if analytic:
                # This expression is known analytically
                new_solved[new_depv] = rel_in_old.subs(inv_trnsfm)
            else:
                # Ok, so it is interconnected to non-analytic dependent variables
                expr_in_old = expr_in_old.subs(eqsd) # df/dt -> x**2+x+...
                expr_in_new = expr_in_old.subs(inv_trnsfm) # put the new variables in place
                new_f[new_depv] = expr_in_new

        new_instance = self.__class__(self, f=new_f, _solved=new_solved)
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
        for depv, (old_expr, sol_symbs) in self._solved.iteritems():
            new_expr = old_expr / sympy.diff(expr_in_old_indep,
                                             self.indepv, order)
            new_expr = new_expr.subs({self.indepv: sympy.solve(
                new_expr - new_indep_symb, self.indepv)[0]})
            new_expr = new_expr.subs({expr_in_old_indep: new_indep_symb})
            new_solved[depv] = new_expr, sol_symbs
        return self.__class__(self, f=new_f, _solved=new_solved,
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

    # Overwrite title for use in e.g. plots
    title = 'System of ordinary differential equations'


    def __init__(self, *args, **kwargs):
        self.param_tokens = self.param_tokens or []
        super(SimpleFirstOrderODESystem, self).__init__(*args, **kwargs)


    def init_f(self):
        # First we need to set the keys (needed when self.expressions()
        # makes look-ups)
        self.f = OrderedDict(
            [(self.mk_func(tok), None) for\
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
            self.param_symbs = [sympy.Symbol(k, real=self.real) for k in\
                                self.param_tokens]
