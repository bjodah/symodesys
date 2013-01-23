
# TODO: Add van der Pool oscillator example and add interfaces
# to both GSL (CythonGSL or own wrapper?) and SUNDIALS (PySundials or own wrapper?)
# TODO: See if it is best to subclass sympy's codegen or use templates with Django template engine.


class IVP_Integrator(object):
    """
    """

    def __init__(self, odesys):
        """

        Arguments:
        - `odesys`:
        - `backend`:
        """
        self._odesys = odesys

    def compile(self):
        """
        To be subclassed.
        Should set self._compiled = True when done
        """
        pass

    def integrate(self, y0, t0, tend, N, abstol = None, reltol = None, h = None):
        if not self._compiled:
            self.compile()
        tout, yout = self.

class GSL_IVP_Integrator()
