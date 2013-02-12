def solve_one_const(f_general_sol, f0, dep_val, const_symb = sympy.symbols('C1')):
    f0_expr = f_general_sol.rhs.subs({dep_val: 0})
    const_val = sympy.solve(sympy.Eq(f0_expr, f0), const_symb)[0]
    return sympy.Eq(f_general_sol.lhs,
                    f_general_sol.rhs.subs({const_symb: const_val}))


class IVP(object):
    """
    Initial Value Problem class

    The class abstracts away the change of initial values
    into parameters when a system of ODEs are reduced analytically
    I.e. the user may still update initial_values, even though in
    ``reality'' it is a parameter which is updated
    """

    def __init__(self, fo_odesys, initial_values, parameters):
        """

        Arguments:
        - `fo_odesys`: First order ODE System
        - `initial_values`: Dictionary mapping dep. var names to initial values
        """
        self._ori_fo_odesys = fo_odesys
        self._fo_odesys = fo_odesys
        self._init_vals = initial_values
        self._init_val_symbs = {}


    def mk_init_val_symb(self, y):
        new_symb = sympy.symbols(str(y) + '_init')
        assert not new_symb == self._fo_odesys.indep_var_symb
        assert not new_symb in self._fo_odesys.dep_var_func_symbs
        assert not new_symb in self._fo_odesys.param_symbs
        return new_symb


    def recursive_analytic_reduction(self):
        """
        Attempts to solve some y's analytically
        """
        solved = {}
        new_init_val_param_symbs = []
        changed_last_loop = True
        while changed_last_loop:
            changed_last_loop = False
            for yi, expr in self._fo_odesys.f.iteritems():
                if yi in solved: continue
                # Check to see if it only depends on itself
                expr = expr.subs(solved)
                Jacobian_row_off_diag = [expr.diff(m) for m \
                                         in self.dep_var_func_symbs if m != yi]
                if all([c == 0 for c in Jacobian_row_off_diag]):
                    # Attempt solution
                    rel = sympy.Eq(yi.diff(x), expr)
                    sol = sympy.dsolve(rel, yi)
                    # ASSIGN NEW SYMBOL TO INITAL VALUE
                    new_symb = (yi, self.mk_init_val_symb(yi))
                    new_init_val_param_symbs[yi] = new_symb
                    sol_init_val = solve_one_const(sol, new_symb, x)
                    solved[yi] = sol_init_val.rhs
                    changed_last_loop = True

        # Post processing of solutions
        self._solved = solved
        for solved_y in solved.keys():
            self._fo_odesys.f.pop(solved_y)
            self._fo_odesys.dep_var_func_symbs.pop(
                self._fo_odesys.dep_var_func_symbs.where(solved_y))

        if len(solved) > 0:
            new_f = {}
            for k, v in self._fo_odesys.f.iteritems():
                new_f[k] = v.subs(solved)
            self._fo_odesys.f = new_f

        self._init_val_symbs.update(new_init_val_param_symbs)
        return new_init_val_param_symbs

    def update_initial_values(self, initial_values):
        """

        """
        for k, v in initial_values.iteritems:
            if k in self._init_vals:
                self._init_vals[k] = v:
            else:
                raise KeyError('Initial value for {} does not exist'.format(k))


    def wrap_integrator(self, Integrator):
        intr = Integrator(self._fo_odesys, self._init_val_symbs)
        def integrate(y0, t0, tend, N, abstol = None, reltol = None, h = None):
            for k in self._init_val_symbs.keys():
                y0.pop(k)
            intr.integrate(y0, t0, tend, N, abstol = None, reltol = None, h = None)
            intr.yout###

        intr.integrate = integrate

