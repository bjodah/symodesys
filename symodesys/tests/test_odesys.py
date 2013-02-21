import sympy

from symodesys.odesys import ODESystem

def test_ODESystem():
    t = sympy.Symbol('t', real = True)
    f = [fnc(t) for fnc in sympy.symbols('f:3', cls = sympy.Function)]
    lmbd = sympy.symbols('lambda:2', real = True)
    k = sympy.Symbol('k', real = True)

    # Coupled decay of 2 variables, the second one controls damping of an oscillator
    odeqs = [sympy.Eq(f[0].diff(t), -lmbd[0] * f[0]),
             sympy.Eq(f[1].diff(t),  lmbd[0] * f[0] - lmbd[1] * f[1]),
             sympy.Eq(f[2].diff(t, 2), -k * f[2] - f[1] * f[2].diff(t))]


    odesys = ODESystem.from_list_of_eqs(odeqs)

    odesys.reduce_to_sys_of_first_order()

    sympy.pprint(odesys.eqs)

if __name__ == '__main__':
    test_ODESystem()
