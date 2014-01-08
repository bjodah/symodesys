# stdlib imports
from collections import namedtuple, OrderedDict, defaultdict

# other imports
import numpy as np
import sympy

from symvarsub.utilities import MaybeRealFunction, reassign_const, get_new_symbs


# project imports
from symodesys.helpers import subs_set, get_without_piecewise


#This might be implemented+named better (unclear if
#  namedtuple appropriate/adds value):
ODEExpr = namedtuple('ODEExpr', 'order expr')

# Abbrevations used in the code:
# ==============================
# depv:   dependent variable (sometimes denoted `y`, as in dydt)
# indepv: independent variable (sometimes denoted `t`, as in dydt)
# f:      dydt denotes the derivates (with their expressions) of
#         dependent variables with respect to the independent variable
# na:     non-analytic
#


class _ODESystemBase(object):
    """
    Central to the _ODESystemBase if the attribute _odeqs
    which should be an OrderedDict(). The keys in the
    dictionary should be of the form:
       sympy.Function('name_of_dependent_variable')(self.indepv)

    User's should use or subclass either FirstOrderODESystem
    or AnyOrderODESystem. _ODESystemBase defines common attributes
    and methods of those derived classes.
    """

    # By default Symbols are taken to represent real valued
    # variables, override this by changing `real` to False:
    # The real attribute makes a difference when solving
    # differential equations using dsolve
    real = True

    # For tracking what dep var func symbs are generated upon
    # order reduction (instantiate as list)
    # it is useful if one does not want to e.g. plot them
    _frst_red_hlprs = None

    # We need holder for analytically solved parts
    # (instantiate to dict):
    _solved = None
    _solved_undefined = None
    _solved_params = None # parameters not used any more

    indepv = None # ODE implies 1 indep. variable,
                  #  set to sympy.Symbol(...)

    param_symbs = None


    #_attrs_to_cmp_for_eq is used for checking equality of instances
    # of _ODESystemBase and its subclasses, the latter requires the
    # list of attrs be a common subset of attributes
    _attrs_to_cmp_for_eq = ['_odeqs', 'param_symbs', '_solved',
                            '_solved_undefined', 'all_depv']

    # _canonical_attrs specifies what attributes are the _defining_
    # attributes of the (sub)class. It maybe changed in subclasses.
    _canonical_attrs = ['_odeqs', 'indepv', 'param_symbs',
                        '_solved', '_solved_undefined', '_frst_red_hlprs']

    # May be passed as kwargs to __init__ but must not be
    _redundant_attrs = []

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
        self._solved = self._solved or OrderedDict()
        # note: _solved must be OrderedDict for solving
        # of constants to be efficient (introducing one unknown at a time)
        self._solved_undefined = self._solved_undefined or []
        self._solved_params = self._solved_params or []
        self._frst_red_hlprs = self._frst_red_hlprs or []


    def __getitem__(self, key):
        if isinstance(key, sympy.Basic):
            match = None
            for known_symb in self.known_symbs:
                if str(known_symb) == str(key):
                    if match == None:
                        match = known_symb
                    else:
                        # This place should never be reached
                        raise RuntimeError(
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

    def mk_depv(self, key):
        """
        Returns an sympy.Function instance with name key and correct
        dependent variable, if just a symbol is wanted: see mk_symb()
        """
        # metaclasses and special behaviour of subclassed
        # classes of sympy.Function makes this tricky.. see tests
        return MaybeRealFunction(key, [self.indepv], real=self.real)
        #return sympy.Function(key)(self.indepv)


    def mk_symb(self, name):
        """
        Convenience function, if a dependent variable is wanted:
        see mk_depv()
        """
        return sympy.Symbol(name, real=self.real)

    def mk_symb_from_depv(self, depv, tail=None):
        if tail == None: tail = ''
        assert depv in self.all_depv
        return self.mk_symb(depv.func.__name__+tail)


    def ensure_dictkeys_as_symbs(self, val_by_token):
        """
        Convenience function for converting dicts with keys of form
        'y1', 'y2' into sympy.Symbol('y1') through __getitem__ (which
        also checks existance of y1, y2... etc.)
        """
        return {self[k]: v for k, v in val_by_token.items()}


    @property
    def all_depv(self):
        """ Returns a list of all dependent variables of the system """
        return self._odeqs.keys()


    @property
    def na_depv(self):
        """
        Returns a list of all dependent variables of the system
        that have not been solved analytically
        """
        return [x for x in self.all_depv if not x in self._solved]


    @property
    def analytic_depv(self):
        """
        Returns a list of all dependent variables of the system
        that have been solved analytically
        """
        return [x for x in self.all_depv if x in self._solved]


    @property
    def analytic_relations(self):
        """
        Convenience attribute, returns all the expressions of the
        analytic solutions corresponding to the dependent variables
        which have been solved analytically.
        """
        return [self._solved[yi][0] for yi in self.analytic_depv]


    @property
    def analytic_sol_symbs(self):
        """
        Returns a list of symbols introduced when solving
        the expressions of the dependent variables analytically.
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
        Convenience attribute, returns a list of all Symbol/Function
        instances in use in the system.
        """
        return [self.indepv] + self.all_depv + self.param_and_sol_symbs

    # @property
    # def forbidden_symbs(self):
    #     """
    #     Extends self.known_symbs to symbols with names coinciding with
    #     e.g. dependent variables.
    #     """
    #     return self.known_symbs + [self.mk_symb(dvfs.func.__name__) for \
    #                                dvfs in self.all_depv]


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
                self.transform_indepv(value, key)
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
                raise RuntimeError("This point should never be reached")


    @property
    def is_first_order(self):
        """
        Is True if the highest order of the ODE's in the system is 1.
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
        Is True if the system is autonomous (the independent variable is abscent
        in all expressions for the derivatives of the dependent variables)
        """
        for depv, (order, expr) in \
                self._odeqs.iteritems():
            if self.unfunc_depv_in_expr(expr).diff(self.indepv) != 0:
                return False
        return True


    def unfunc_depv_in_expr(self, expr):
        """
        Convenience method, transforms the Function instances of the dependent variables
        possibly present in provided argument `expr` into Symbol instances.
        """
        unfunc_subs = {dvfs: self.mk_symb(dvfs.func.__name__) for \
                       dvfs in self.all_depv}
        return expr.subs(unfunc_subs)


    def refunc_depv_in_expr(self, expr):
        """ The inverse of unfunc_depv_in_expr """
        refunc_subs = {self.mk_symb(dvfs.func.__name__): dvfs for \
                       dvfs in self.all_depv}
        return expr.subs(refunc_subs)


    @property
    def is_homogeneous(self):
        """ TODO implement this """
        raise NotImplementedError


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
        Returns a sympy.Eq for diff eq of ``depv'',
        if it is solved analytically, the analytic solution
        is returned
        """
        order, expr = self._odeqs[depv]
        if depv in self._solved:
            return sympy.Eq(depv, self._solved[depv][0])
        else:
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
                if not all([wrt[0] == x for x in wrt]):
                    raise NotImplementedError('Systems of partial'+\
                            ' differential equations not supported!')
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


    def attempt_analytic_sol(self, depv, hypoexpr, sol_consts=None):
        """
        Checks if provided analytic (similiar to
        sympy.solvers.ode.checkodesol) expr ``hypoexpr'' solves
        diff eq of ``depv'' in odesys. If it does new symbols are
        generated and saved toghether with hypoexpr in self._solved
        """
        if depv in self._solved:
            raise ValueError("{} already solved analytically".format(dpev))
        else:
            eq = self.eq(depv)

        if bool(eq.doit()):
            if sol_consts:
                assert len(get_new_symbs(
                    hypoexpr, self.known_symbs + sol_consts)) == 0
            else:
                sol_consts = get_new_symbs(
                    hypoexpr, self.known_symbs + sol_consts)
            self._solved[depv] = hypoexpr, sol_consts

            # Now express the new symbols in terms of the new ones
            if order > 1:
                eq = self.eq(depv).subs({depv: hypoexpr})

            eq0 = eq.subs({self.indepv: t0})
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
            (k, ODEExpr(1, v.expr)) for k, v in self._odeqs.items()])

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
            candidate = self.mk_depv(fnc.func.__name__ + '_h' + str(i))
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
        _frst_red_hlprs = self._frst_red_hlprs[:]
        if default_red_init_val == None: default_red_init_val = 0.0

        # TODO, revise _frst_red_hlprs format (cur.  len 3 tuple)
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
            _frst_red_hlprs.extend([
                (fnc, o, expr) for o, expr in hlpr.iteritems()])
        new_instance = FirstOrderODESystem(
            f=OrderedDict([(k, v.expr)for k,v in new_odeqs.items()]),
            indepv=self.indepv, param_symbs=self.param_symbs[:])
        new_instance._frst_red_hlprs = _frst_red_hlprs
        if y0:
            for source, order, helper in _frst_red_hlprs:
                if not helper in y0:
                    y0[helper] = default_red_init_val
        return new_instance


class FirstOrderODESystem(_ODESystemBase):
    """
    This class provides an efficient interface to the manipulation and
    export of first order ODE Systems, it is a special case in the
    sense that it only contains first order derivatives (which is a
    requirement for the numerical treatment of ODE Systems).

    The attribute `f` is given a special role. It is a simpliciation
    of _odeqs attribute which now has a redundant order part of its
    (order, expr) tuples.  Furthermore Hence, `f` is a dictionary
    mapping the Symbol of the dependent variable to the expression of
    the first derivative of same dependent variable with respect to
    the independnt variable.

    When the ODE systems equations is generated from user data this
    class is to be subclassed to provide similar routines of the
    SimpleFirstOrderODESystem which is to be used when explicitly
    defining an ODE System.
    """

    f = None

    _canonical_attrs = ['f', 'indepv', 'param_symbs',
                        '_solved', '_solved_undefined', '_frst_red_hlprs']

    _redundant_attrs = ['_odeqs']

    info = None

    def __init__(self, odesys = None, **kwargs):
        # Using property functions makes overwriting harder therefore
        # 'indepv', 'all_depv', 'param_symbs', 'f'
        # are initialized via instance methods that by default
        # initializes empty list / dict only if attribute is missing
        super(FirstOrderODESystem, self).__init__(odesys, **kwargs)
        self._init_param_symbs()
        if self.f == None:
            self._init_f()
        assert self.is_first_order
        self.info = self.info or defaultdict(int)


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

            if self.unfunc_depv_in_expr(expr).diff(self.indepv) != 0:
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
            (k, ODEExpr(1, self.f[k])) for k in self.all_depv])

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


    def recursive_analytic_auto_sol(self, complexity=0, cb=None, logger=None):
        """
        Solves equations one by one (hence it can only find solutions
        for independent equations)

        Input arguments:
        -`complexity`: Positivt integer denoting the degree of complexity
            (an arbitrary measure) to solve for:
          0: solve only constants
          1: try to solve only differentials in independent variable
          2: try to solve only differentials in auto-dependent variable
          3: try to solve only differentials in independent and/or auto-dependent
             variables
        -`cb`: a callback with signature (expr, depv, const_symb) for
               custom assignment of integration constants. The callback
               returns a new expr

        TODO: sometimes the different solutions are valid for
        different conditions of choosen parameters.
        In sympy it indicated by the solution consisting of multiple
        sympy.Piecewise instances. See docstring of
        FirstOrderODESystem._dsolve for more information.
        """
        nsolved0 = len(self._solved)
        changed_last_loop = True
        while changed_last_loop:
            changed_last_loop = False
            for yi, expr in self.f.iteritems():
                if yi in self._solved: continue
                # Substitute solved differentials of dependent variables
                expr = expr.subs({yi: expr for yi, (expr, sol_symbs) \
                                  in self._solved.iteritems()})
                if complexity == 0:
                    Jac_row = [expr.diff(m) for m \
                               in self.na_depv]
                    if all([c == 0 for c in Jac_row]):
                        # no depv in expr
                        if not expr.has(self.indepv):
                            # It must be a constant!
                            self._dsolve(yi, expr, cb)
                            if logger:
                                logger.info('Solved {}={}'.format(
                                    yi, self._solved[yi][0]))
                            changed_last_loop = True

                if complexity == 1 or complexity == 3:
                    Jac_row = [expr.diff(m) for m \
                               in self.na_depv]
                    if all([c == 0 for c in Jac_row]):
                        # no depv in expr
                        self._dsolve(yi, expr, cb)
                        if logger:
                            logger.info('Solved {}={}'.format(
                                yi, self._solved[yi][0]))
                        changed_last_loop = True

                if complexity == 2 or complexity == 3:
                    # Check to see if it only depends on itself
                    Jac_row_off_diag = [expr.diff(m) for m \
                                        in self.na_depv if m != yi]
                    if all([c == 0 for c in Jac_row_off_diag]):
                        # expr is auto-dependent
                        self._dsolve(yi, expr, cb)
                        if logger:
                            logger.info('Solved auto-dependent {}={}'.format(
                                yi, self._solved[yi][0]))
                        changed_last_loop = True
        return len(self._solved) - nsolved0


    def _dsolve(self, depv, expr, cb=None):
        """
        Run sympy's dsolve routine on the expression

        Derivative(`depv`, self.indepv) = `expr`

        provide callback `cb` for giving the integration constant
        a custom expression. cb has the signature:
        cb(expr, depv, const_symb) and returns a tuple:
        new_expr, new_const_symb

        Currently if sympy.dsolve returns a PieceWise result,
        only the default expression is stored while the other conditions are
        stored as list of illegal conditions. This is done not to have
        a branching tree of analytic solutions. If the interest is in e.g.:
        param_a == param_b then consider ODESys.subs({param_b: param_a}) as
        a possible work around (with the obvious caveat that the ODE System
        is significantly less flexible from that point on).
        """
        # Attempt solution (actually: assume success)
        rel = sympy.Eq(depv.diff(self.indepv), expr)
        sol = sympy.dsolve(rel, depv)

        # If sol contains a Piecewise definition,
        # accept the default solution and store
        # the others as undefined cases.
        if sol.rhs.has(sympy.Piecewise):
            sol_expr, undefined_cases = get_without_piecewise(sol.rhs)
            self._solved_undefined.extend(undefined_cases)
        else:
            sol_expr = sol.rhs

        # Assign new symbol to inital value
        sol_expr, rea, not_rea = reassign_const(
            sol_expr, depv.func.__name__+'C', self.known_symbs)
        assert len(not_rea) == 0

        if cb:
            const_symb_in_init, new_const = cb(sol_expr, depv, rea[0])
            assert len(rea) == 1
            sol_expr = sol_expr.subs({
                rea[0]: const_symb_in_init})
            self._new_solve(depv, sol_expr, [new_const], expr)
        else:
            self._new_solve(depv, sol_expr, rea, expr)


    def _new_solve(self, depv, sol_expr, sol_symbs, ori_expr):
        """
        Manipulates self._solved and handles param_symbs
        """
        # Handle param_symbs
        ## Identify what params are used in current expr
        exclusively_present = []
        for symb in self.param_symbs:
            if symb in ori_expr:
                for dv in self.all_depv:
                    if dv != depv:
                        if dv in self._solved:
                            if symb in self._solved[dv][0]:
                                break
                        else:
                            if symb in self._odeqs[dv][1]:
                                break
                else: # exlusive
                    exclusively_present.append(symb)
        assert len(sol_symbs) >= len(exclusively_present)
        self._solved_params.extend(exclusively_present)
        self._solved[depv] = sol_expr, sol_symbs


    def _init_param_symbs(self):
        """
        To be subclassed (or add attribute (list): param_symbs)

        should set self.param_symbs to a list of
        sympy.symbols(``token_string'') instances.
        The order in this list defines indices in vectors
        and matrices used by underlying numerical integration.  (When
        subclassing, sympy.symarray might be useful.)
        """
        if self.param_symbs == None:
            self.param_symbs = []


    def _init_f(self):
        """
        To be subclassed.

        *Is only exectuted if and only if self.f != None
        *self._init_f() must:
          set self.f to a OrderedDict with the first-order
          derivatives as values (and dependent variable sympy.Function
          instances as keys, use mk_depv to create)
        """
        self.f = OrderedDict()


    @property
    def analytic_subs(self):
        return {depv: expr for depv, (expr, symbs) \
                in self._solved.items()}


    @property
    def na_f(self):
        if len(self._solved) > 0:
            return OrderedDict(
                [(k,  self.f[k].subs(self.analytic_subs)) for k in \
                 self.na_depv])
        else:
            return self.f


    @property
    def jac(self):
        fmat = sympy.Matrix(
            1, len(self.all_depv),
            lambda q, i: self.f[self.all_depv[i]]
            )
        return fmat.jacobian(self.all_depv)


    @property
    def na_jac(self):
        na_fmat = sympy.Matrix(
            1, len(self.na_depv),
            lambda q, i: self.f[self.na_depv[i]]
            )
        return na_fmat.jacobian(self.na_depv).subs(self.analytic_subs)


    def _get_num_subs(self, depv_vals, param_vals):
        depv_subs = dict(zip(self.na_depv, depv_vals))
        param_subs = dict(zip(self.param_symbs, param_vals))
        subs = {}
        for d in [depv_subs, param_subs]:
            subs.update(d)
        return subs


    def param_val_lst(self, param_vals_by_symb):
        return [param_vals_by_symb[k] for k in self.param_symbs]


    def _subs_depv_dict(self, d, indepv_val, depv_vals_d, param_vals_d):
        """
        Provide e.g. self.f or self.na_f as `d`.
        See self.subs_f
        """
        return OrderedDict([
            (depv, expr.subs(depv_vals_d).subs(param_vals_d).subs(
                {self.indepv: indepv_val})) for depv, expr in \
            d.items()])


    def subs_f(self, indepv_val, depv_d, param_d):
        """
        Returns dydt (f) with substituted expressions using the input
        values:

        -`indepv_val`: scalar
        -`depv_d`: dictionary
        -`param_d`: dictionary
        """
        return self._subs_depv_dict(self.f, indepv_val, depv_d, param_d)


    def subs_na_f(self, indepv_val, depv_d, param_d):
        """ Same as self.subs_f, but only for the dependent variables
        not solved analytically, and param_d now also include
        integration constants from the analytic treatment """
        return self._subs_depv_dict(self.na_f, indepv_val,
                               depv_d, param_d)


    def evaluate_f(self, indepv_val, depv_arr, param_arr):
        """
        Convenience function for evaluating dydt (f) for
        provided input values:
        -`indepv_val`: scalar
        -`depv_arr`: array of floats (ordered as self.all_depv)
        -`param_arr`: array of floats (ordered as self.param_symbs)

        It returs an array with the numerical value of the derivatives
        ordered as self.all_depv
        """
        assert len(depv_arr) == len(self.all_depv)
        assert len(param_arr) == len(self.param_symbs)
        return self.subs_f(indepv_val,
                    dict(zip(self.all_depv, depv_arr)),
                    dict(zip(self.param_symbs, param_arr))).values()


    def evaluate_na_f(self, indepv_val, depv_arr, param_arr):
        """
        Same as evaluate_f but only for dependent vairables whose
        derivates have not been integrated analytically, and param_arr
        should now also include integration constants.
        """
        assert len(depv_arr) == len(self.na_depv)
        assert len(param_arr) == len(self.param_and_sol_symbs)
        self.info['nfev'] += 1 # SciPy integrator callback
        return self.subs_na_f(
            indepv_val,
            dict(zip(self.na_depv, depv_arr)),
            dict(zip(self.param_and_sol_symbs, param_arr))).values()


    def _subs_list_of_lists(self, ll, indepv_val, depv_d, param_d):
        return [[cell.subs(param_d).subs(depv_d).subs(
            {self.indepv: indepv_val}) for \
                 cell in row] for row in ll]


    def subs_jac(self, indepv_val, depv_d, param_d):
        """
        Returns the Jacobian dfdy (f === dydt) with substituted
        expressions using the input values:

        -`indepv_val`: scalar
        -`depv_d`: dictionary
        -`param_d`: dictionary
        """
        return self._subs_list_of_lists(self.jac.tolist(),
                                        indepv_val, depv_d, param_d)


    def subs_na_jac(self, indepv_val, depv_d, param_d):
        return self._subs_list_of_lists(self.na_jac.tolist(),
                                        indepv_val, depv_d, param_d)


    def evaluate_jac(self, indepv_val, depv_arr, param_arr):
        """
        Convenience function for evaluating the Jacobian dfdy
        (f === dydt) for provided input values:
        -`indepv_val`: scalar
        -`depv_arr`: array of floats (ordered as self.all_depv)
        -`param_arr`: array of floats (ordered as self.param_symbs)

        It returs an rank 2 array with the numerical value of the
        Jacobian elements ordered as self.all_depv
        """
        assert len(depv_arr) == len(self.all_depv)
        assert len(param_arr) == len(self.param_symbs)
        return self.subs_jac(
            indepv_val,
            dict(zip(self.all_depv, depv_arr)),
            dict(zip(self.param_symbs, param_arr)))



    def evaluate_na_jac(self, indepv_val, depv_arr, param_arr):
        assert len(depv_arr) == len(self.na_depv)
        assert len(param_arr) == len(self.param_and_sol_symbs)
        self.info['njev'] += 1 # SciPy integrator callback
        return self.subs_na_jac(
            indepv_val,
            dict(zip(self.na_depv, depv_arr)),
            dict(zip(self.param_and_sol_symbs, param_arr)))


    @property
    def dfdt(self):
        f_subs = {k.diff(self.indepv): v for k, v in \
             self.f.iteritems()}
        return OrderedDict([
            (y, self.refunc_depv_in_expr(
                self.unfunc_depv_in_expr(self.f[y]).diff(
                    self.indepv).subs(f_subs))) for \
            y in self.all_depv])


    @property
    def na_dfdt(self):
        na_f_subs = {k.diff(self.indepv): v for k, v in \
                     self.na_f.iteritems()}
        return OrderedDict([
            (y, self.refunc_depv_in_expr(
                self.unfunc_depv_in_expr(self.na_f[y]).diff(
                    self.indepv).subs(na_f_subs))) for \
            y in self.na_depv])


    def subs_dfdt(self, indepv_val, depv_d, param_d):
        """
        Returns dfdt (d2ydt2) with substituted expressions using the input
        values:

        -`indepv_val`: scalar
        -`depv_d`: dictionary
        -`param_d`: dictionary
        """
        return self._subs_depv_dict(self.dfdt, indepv_val, depv_d, param_d)


    def subs_na_dfdt(self, indepv_val, depv_d, param_d):
        """ Same as self.subs_dfdt, but only for the dependent variables
        not solved analytically, and param_d now also include
        integration constants from the analytic treatment """
        return self._subs_depv_dict(self.na_dfdt, indepv_val,
                                    depv_d, param_d)


    def evaluate_dfdt(self, indepv_val, depv_arr, param_arr):
        """
        Convenience function for evaluating dfdt (d2ydt2) for
        provided input values:
        -`indepv_val`: scalar
        -`depv_arr`: array of floats (ordered as self.all_depv)
        -`param_arr`: array of floats (ordered as self.param_symbs)

        It returs an array with the numerical value of the derivatives
        ordered as self.all_depv
        """
        assert len(depv_arr) == len(self.all_depv)
        assert len(param_arr) == len(self.param_symbs)
        return self.subs_dfdt(
            indepv_val,
            dict(zip(self.all_depv, depv_arr)),
            dict(zip(self.param_symbs, param_arr))).values()


    def evaluate_na_dfdt(self, indepv_val, depv_arr, param_arr):
        """
        Same as evaluate_dfdt but only for dependent vairables whose
        derivates have not been integrated analytically, and param_arr
        should now also include integration constants.
        """
        assert len(depv_arr) == len(self.na_depv)
        assert len(param_arr) == len(self.param_and_sol_symbs)
        self.info['ndfdtev'] += 1 # SciPy integrator callback
        return self.subs_na_dfdt(
            indepv_val,
            dict(zip(self.na_depv, depv_arr)),
            dict(zip(self.param_and_sol_symbs, param_arr))).values()


    def evaluate_d2ydt2(self, indepv_val, depv_arr, param_arr):
        return np.dot(
            self.evaluate_jac(indepv_val, depv_arr, param_arr),
            self.evaluate_f(indepv_val, depv_arr, param_arr)
        ) + self.evaluate_dfdt(indepv_val, depv_arr, param_arr)


    def evaluate_na_d2ydt2(self, indepv_val, depv_arr, param_arr):
        return np.dot(
            self.evaluate_na_jac(indepv_val, depv_arr, param_arr),
            self.evaluate_na_f(indepv_val, depv_arr, param_arr)
        ) + self.evaluate_na_dfdt(indepv_val, depv_arr, param_arr)


    def is_stiff(self, indepv_val, depv_vals, param_vals,
                 indepv_val_end, criteria=100):
        """
        TODO: Add reference to definition of stiffness
        """
        num_jac = self.dydt_jac(indepv_val, depv_vals, param_vals)
        eigen_values = scipy.linalg.eigvals(num_jac)
        mineigen = np.min(eigenvalues)
        tspan = tend-t0
        if mineigen < 0:
            if -mineigen*tspan > criteria:
                return True
        return False


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

        # d/dt(old): expr_in_old
        eqsd = {eq.lhs: eq.rhs for eq in self.eqs}
        new_solved = {}
        for new_depv, rel_in_old in trnsfm.iteritems():
            # First use our expressions for solved variables
            rel_in_old = rel_in_old.subs(dict(zip(
                self.analytic_depv, self.analytic_relations)))
            # what does the derivatives look like?
            expr_in_old = rel_in_old.diff(self.indepv)
            analytic = True
            for depv in self.na_depv:
                if expr_in_old.has(depv.diff(self.indepv)):
                    analytic = False
                    break
            if analytic:
                # This expression is known analytically
                new_solved[new_depv] = rel_in_old.subs(inv_trnsfm)
            else:
                # Ok, so it is interconnected to non-analytic ,
                # dependent variables
                # df/dt -> x**2+x+...
                expr_in_old = expr_in_old.subs(eqsd)

                # put the new variables in place
                expr_in_new = expr_in_old.subs(inv_trnsfm)
                new_f[new_depv] = expr_in_new

        new_instance = self.__class__(self, f=new_f,
                                      _solved=new_solved)
        return new_instance


    def transform_indepv(self, new_indepv_symb, expr_in_old_indepv):
        new_f = OrderedDict()
        for depv, old_expr in self.f.iteritems():
            new_expr = old_expr / \
                       sympy.diff(expr_in_old_indepv, self.indepv)
            new_expr = new_expr.subs({self.indepv: sympy.solve(
                new_expr - new_indepv_symb, self.indepv)[0]})
            new_expr = new_expr.subs({expr_in_old_indepv: new_indepv_symb})
            new_f[depv] = new_expr

        # Handle solved:
        new_solved = {}
        for depv, (old_expr, sol_symbs) in self._solved.iteritems():
            new_expr = old_expr / sympy.diff(expr_in_old_indepv,
                                             self.indepv, order)
            new_expr = new_expr.subs({self.indepv: sympy.solve(
                new_expr - new_indepv_symb, self.indepv)[0]})
            new_expr = new_expr.subs({expr_in_old_indepv: new_indepv_symb})
            new_solved[depv] = new_expr, sol_symbs
        return self.__class__(self, f=new_f, _solved=new_solved,
                              indepv=new_indepv_symb)


class SimpleFirstOrderODESystem(FirstOrderODESystem):
    """
    This class provides convenience methods for generating the
    symbols of the idependent variable symbol, dependent variable symbols
    and parameter symbols. It is useful when the equations are not
    algorithmatically generated but by user subclassing (of this class).

    Essentially, ``tokens'' are then the generating strings of the
    symbol instances of parent class.
    """

    indepv = sympy.symbols('t') # ODE implies 1 independent variable

    # The following must be provided in this "Simple" subclass
    # (string reprs instead of sympy.symbols)
    depv_tokens = None
    param_tokens = None
    expressions = None

    # Overwrite title for use in e.g. plots
    title = None


    def __init__(self, *args, **kwargs):
        self.param_tokens = self.param_tokens or []
        super(SimpleFirstOrderODESystem, self).__init__(*args, **kwargs)


    def _init_f(self):
        # First we need to set the keys (needed when self.expressions()
        # makes look-ups)
        self.f = OrderedDict(
            [(self.mk_depv(tok), None) for\
                              tok in self.depv_tokens])
        for tok in self.depv_tokens:
            self.f[self[tok]] = self.expressions[self[tok]]


    def param_vals_by_symb(self, params_by_token):
        return dict([(self[k], params_by_token[k]) for\
                     k in self.param_tokens])


    def _init_param_symbs(self):
        if self.param_symbs == None:
            self.param_symbs = [self.mk_symb(k) for k \
                                in self.param_tokens]
