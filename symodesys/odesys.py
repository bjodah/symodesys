
# Experimental
import sympy

class ODESystem(object):

    _1st_ordr_red_helper_fncs = []

    def __init__(self, odes, indep_var_symb, parameters):
        self._odes = odes
        self._indep_var_symb = indep_var_symb
        self._parameters = parameters

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
            while candidate in self._odes:
                candidate = candidate + '_h'
            helpers[i] = candidate
        return helpers

    def reduce_to_sys_of_first_order(self):
        for fnc, (order, expr) in self._odes.iteritems():
            if order == 1: continue
            hlpr = self.get_helper_fnc(fnc, order)
            self._odes[fnc] = (1, hlpr[1])
            for o in range(1, order - 1):
                self._odes[hlpr[o]] = (1, hlpr[o + 1])
            subsd = {sympy.Derivative(fnc, self._indep_var_symb, i): hlpr[i] \
                          for i in range(1, order)}
            self._odes[hlpr[order - 1]] = (1, expr.subs(subsd))
            self._1st_ordr_red_helper_fncs.extend([(fnc, o, expr) for o, expr \
                                                   in hlpr.iteritems()])
            # self._odes changed > call this function recursively
            self.reduce_to_sys_of_first_order()
            break

    @classmethod
    def from_list_of_eqs(cls, lst):
        """
        Determine independnt variable, dependent variables
        and parameter variables and perform Sanity checks
        """
        fncs = []
        odes = {} # independet variable as key, (order, expr) as value
        indep_var_symb = None
        parameters = set()
        for eq in lst:
            assert eq.lhs.is_Derivative
            fnc, wrt = eq.lhs.args[0], eq.lhs.args[1:]
            assert fnc not in odes
            if indep_var_symb == None:
                assert all([wrt[0] == x for x in wrt]) # No PDEs!
                indep_var_symb = wrt[0]
            else:
                assert all([indep_var_symb == x for x in wrt])
            odes[fnc] = (len(wrt), eq.rhs)

        for fnc, (order, expr) in odes.iteritems():
            # Sanity check
            for check_fnc in odes.keys():
                # Make sure expr don't contain derivatives larger than in odes order:
                if expr.has(check_fnc):
                    for arg in expr.args:
                        if arg.has(check_fnc):
                            if arg.is_Derivative:
                                fnc, wrt = args[0], args[1:]
                                assert not len(wrt) > odes[check_fnc][0]

            for atom in expr.atoms():
                if not atom in odes.keys() and atom != indep_var_symb and not atom.is_Number:
                    parameters.add(atom)

        return cls(odes, indep_var_symb, parameters)

    @property
    def eqs(self):
        """
        Returns a list of Sympy Eq instances describing the ODE system
        """
        return [sympy.Eq(k.diff(self._indep_var_symb, v[0]), v[1]) for k, v in self._odes.iteritems()]
