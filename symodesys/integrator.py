
# TODO: Add van der Pool oscillator example and add interfaces
# to both GSL (CythonGSL or own wrapper?) and SUNDIALS (PySundials or own wrapper?)
# TODO: See if it is best to subclass sympy's codegen or use templates with Django template engine.


class Integrator(object):
    """

    TODO: Check if it is valuable to have taylor expanded version of approximated d/dC (Jac)
    TODO: Check if it is valuable to have (taylor expanded version of approximated) inv(Jac)
    """

    def __init__(self, reaction_system, C0, abstol = None, reltol = None, h = None, logarithmic = False):
        """

        Arguments:
        - `reaction_system`:
        - `C0`:
        - `h`:
        - `abstol`:
        - `reltol`:
        - `logarithmic`:
        """
        self._reaction_system = reaction_system
        self._C0 = C0
        self._abstol = abstol
        self._reltol = reltol
        self._h = h
        self._logarithmic = logarithmic

    def write_code(self, dest):
        pass

    def compile_code(self):
        pass

    def integrate(self, t0, tend):
        pass


