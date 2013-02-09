

class IVP(Problem):
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
        self._initial_values = initial_values
        self.recursive_analytic_reduction()

    def attempt_analytic_reduction(self, fo_odesys):
        """
        Attempt to solve part of first order ode sys
        analytically
        """
        pass


    def recursive_analytic_reduction(self):
        fo_odesys = self._fo_odesys
        old_fo_odesys = None
        while (old_fo_odesys != fo_odesys):
            old_fo_odesys = fo_odesys
            fo_odesys, = self.attempt_analytic_reduction(fo_odesys)
        if self._fo_odesys != fo_odesys:
            self._fo_odesys = fo_odesys


    def update_initial_values(self, initial_values):
        """

        """
        for k, v in initial_values.iteritems:
            if k in self._initial_values:
                self._initial_values[k] = v:
            else:
                raise KeyError('Initial value for {} does not exist')

