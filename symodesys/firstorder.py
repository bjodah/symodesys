
class FirstOrderODESystem(object):
    """
    When the ODE systems equations is generated from user
    data this class is to be subclassed to provide the routines
    described below
    """
    indep_var_symb = sympy.symbols('t') # ODE implies 1 indep. variable

    ny = None # Size of ODE system
    nk = None # Number of parameters

    _dep_var_basesymb = 'y'
    _param_basesymb = 'k'

    @property
    def dep_var_symbs(self):
        """
        To be subclassed
        should return e.g.  sympy.symarray('y', 3)
        """
        return sympy.symarray(self._dep_var_basesymb, self.ny)
    y = dep_var_symbs

    @property
    def param_symbs(self):
        """
        To be subclassed
        should return e.g.  sympy.symarray('k', 7)
        """
        return sympy.symarray(self._param_basesymb, self.nk)
    k = param_symbs

    @property
    def param_default_values(self):
        """
        To be subclassed
        should return list of same length as self.param_symbs()
        with default values. None represents no default value.
        """
        return np.zeros(self.nk)


    @property
    def f(self):
        """
        To be subclassed, should return a sympy matrix of dimension
        (1, len(self.dep_var_symb))
        for the first-order derivatives of the self.dep_var_symb
        expressed solely in constants, indep_var_symb, dep_var_symb
        and param_symbs
        """
        pass

    def _fmat(self):
        pass

    def transform_indep_var_to_log_scale(self):
        pass

    def transform_dep_vars_to_log_scale(self):
        pass

    def reduce_sys_by_solving_decoupled_vars_analytically(self):
        pass
