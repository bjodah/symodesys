
# Experimental
import sympy

from collections import OrderedDict

class ODESystem(object):

    indep_var_symb = None # ODE implies 1 indep. variable, set to sympy.symbol(...)
    dep_var_func_symbs = None # Set to a list of kind:
                              #    sympy.Function('name')(self.indep_var_symb)
    param_symbs = None
    f = None
    param_vals_by_symb = None

    @property
    def is_autonomous(self):
        for dep_var_func_symb, (order, expr) in self._odeqs_by_indep_var.iteritems():
            unfunc_subs = {dvfs: sympy.symbol(dvfs.func.__name__) for dvfs \
                           in self.dep_var_func_symbs}
            if expr.subs(unfunc_subs).diff(self.indep_var_symb) != 0: return False
        return True

    @property
    def is_linear(self):
        # Most easily done for first order system?
        for dep_var_func_symb, (order, expr) in self._odeqs_by_indep_var.iteritems():
            for wrt in self.dep_var_func_symbs:
                expr = expr.diff(wrt)

            unfunc_subs = {dvfs: sympy.symbol(dvfs.func.__name__) for dvfs \
                           in self.dep_var_func_symbs}
            if expr.subs(unfunc_subs).diff(self.indep_var_symb) != 0: return False

        pass

    @property
    def is_homogeneous(self):
        pass

    @property
    def eqs(self):
        """
        Returns a list of Sympy Eq instances describing the ODE system
        """
        return [sympy.Eq(k.diff(self._indep_var_symb, v[0]), v[1]) for k, v in \
                self._odeqs_by_indep_var.iteritems()]

    def do_sanity_check_of_odeqs(self):
        for fnc, (order, expr) in odeqs_by_indep_var.iteritems():
            # Sanity check
            for check_fnc in odeqs_by_indep_var.keys():
                # Make sure expr don't contain derivatives larger
                # than in odeqs_by_indep_var order:
                if expr.has(check_fnc):
                    for arg in expr.args:
                        if arg.has(check_fnc):
                            if arg.is_Derivative:
                                fnc, wrt = args[0], args[1:]
                                assert not len(wrt) > odeqs_by_indep_var[check_fnc][0]

    def __eq__(self, other):
        for attr in self._attrs_to_cmp_for_eq:
            if getattr(self, attr) != getattr(other, attr): return False
        return True



class AnyOrderODESystem(ODESystem):

    #_attrs_to_cmp_for_eq is used for checking equality of class instances
    _attrs_to_cmp_for_eq = ['_odeqs_by_indep_var', 'indep_var_symb', 'param_symbs']

    # When reducing the order of the system to 1st order
    # a lot of helper variables are introduced. If one do
    # not want to e.g. plot them, the following list if useful:
    _1st_ordr_red_helper_fncs = []


    def __init__(self, odeqs_by_indep_var, indep_var_symb, param_symbs):
        self._odeqs_by_indep_var = odeqs_by_indep_var
        self._indep_var_symb = indep_var_symb
        self._param_symbs = param_symbs

    def attempt_analytic_reduction(self):
        pass
        #return {dep_var: expr}

    def get_helper_fnc(self, fnc, order):
        """
        Returns a list of of length order - 1
        for use in reformulation of higher order
        ODE eq in first order ODE (sub) sys.
        """
        helpers = {}
        for i in range(1, order):
            candidate = sympy.Function(str(fnc.func.__name__) + '_h' + str(i))(self._indep_var_symb)
            while candidate in self._odeqs_by_indep_var:
                candidate = candidate + '_h'
            helpers[i] = candidate
        return helpers

    def reduce_to_sys_of_first_order(self):
        for fnc, (order, expr) in self._odeqs_by_indep_var.iteritems():
            if order == 1: continue
            hlpr = self.get_helper_fnc(fnc, order)
            self._odeqs_by_indep_var[fnc] = (1, hlpr[1])
            for o in range(1, order - 1):
                self._odeqs_by_indep_var[hlpr[o]] = (1, hlpr[o + 1])
            subsd = {sympy.Derivative(fnc, self._indep_var_symb, i): hlpr[i] \
                          for i in range(1, order)}
            self._odeqs_by_indep_var[hlpr[order - 1]] = (1, expr.subs(subsd))
            self._1st_ordr_red_helper_fncs.extend([(fnc, o, expr) for o, expr \
                                                   in hlpr.iteritems()])
            # self._odeqs_by_indep_var changed > call this function recursively
            self.reduce_to_sys_of_first_order()
            break

    @classmethod
    def from_list_of_eqs(cls, lst):
        """
        Determine independnt variable, dependent variables
        and parameter variables and perform Sanity checks
        """
        fncs = []

        # independet variable as key, (order, expr) as value
        odeqs_by_indep_var = OrderedDict()
        indep_var_symb = None
        for eq in lst:
            assert eq.lhs.is_Derivative
            fnc, wrt = eq.lhs.args[0], eq.lhs.args[1:]
            assert fnc not in odeqs_by_indep_var
            if indep_var_symb == None:
                assert all([wrt[0] == x for x in wrt]) # No PDEs!
                indep_var_symb = wrt[0]
            else:
                assert all([indep_var_symb == x for x in wrt])
            odeqs_by_indep_var[fnc] = (len(wrt), eq.rhs)

        param_symbs = set()
        for fnc, (order, expr) in odeqs_by_indep_var.iteritems():
            for atom in expr.atoms():
                if not atom in odeqs_by_indep_var.keys() and \
                       atom != indep_var_symb and not atom.is_Number:
                    param_symbs.add(atom)
        new_instance = cls(odeqs_by_indep_var, indep_var_symb, param_symbs)
        new_instance.do_sanity_check_of_odeqs()
        return new_instance

    @property
    def dep_var_func_symbs(self):
        return self._odeqs_by_indep_var.keys()
