
# Experimental


class ODESystem(object):

    self._1st_ordr_red_helper_fncs = []

    def __init__(self, odes, indep_var, parameters):
        self._odes = odes
        self._indep_var = indep_var
        self._parameters = parameters

    def attempt_analytic_reduction(self):
        pass
        #return {dep_var: expr}


    def reduce_to_sys_of_first_order(self):
        for fnc, (order, expr) in self._odes.iteritems():
            if order == 1: continue
            hlpr_fnc = self.get_helper_fnc(fnc, order)
            self._odes[fnc] = order - 1, hlpr_fnc
            self._odes[hlpr_fnc] = expr.subs(fnc.der(self._indep_var, order - 1), fnc)
            # FIX THIS WITH PEN AND PAPER!!!
            raise NotImplementedError
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
        indep_var = None
        parameters = set()
        for eq in lst:
            assert eq.lhs.is_Derivative
            fnc, wrt = eq.lhs.args[0], eq.lhs.args[1:]
            assert fnc not in odes
            if indep_var == None:
                assert all([wrt[0] == x for x in wrt]) # No PDEs!
                indep_var = wrt[0]
            else:
                assert all([indep_var == x for x in wrt])
            odes[fnc] = (len(wrt), eq.lhs)

        for fnc, (order, expr) in odes.iteritems():
            # Sanity check
            for check_fnc in odes.keys():
                assert not expr.has(check_fnc)

            for atom in expr.atoms():
                if not atom in odes.keys() and atom != indep_var and not atom.is_Number:
                    parameters.add(atom)

        return cls(odes, indep_var, parameters)
