import sympy

from symodesys.odesys import ODESystem, FirstOrderODESystem, SimpleFirstOrderODESystem, AnyOrderODESystem

def test_AnyOrderODESystem___1():
    t = sympy.Symbol('t', real = True)
    f = [fnc(t) for fnc in sympy.symbols('f:3', cls = sympy.Function)]
    lmbd = sympy.symbols('lambda:2', real = True)
    k = sympy.Symbol('k', real = True)

    # Coupled decay of 2 variables, the second one controls damping of an oscillator
    odeqs = [sympy.Eq(f[0].diff(t), -lmbd[0] * f[0]),
             sympy.Eq(f[1].diff(t),  lmbd[0] * f[0] - lmbd[1] * f[1]),
             sympy.Eq(f[2].diff(t, 2), -k * f[2] - f[1] * f[2].diff(t))]


    odesys = AnyOrderODESystem.from_list_of_eqs(odeqs)

    odesys.reduce_to_sys_of_first_order()

    sympy.pprint(odesys.eqs)

def test_ODESystem___2():
    class Decay(SimpleFirstOrderODESystem):
        dep_var_tokens = 'u',
        param_tokens   = 'lambda_u',

        @property
        def expressions(self):
            return {self['u']: self['lambda_u'] * -self['u']}

    d=Decay()

    s = d.mk_func('s')
    #print s._eval_is_real()
    assert s.is_real
    assert d.mk_func('s') == s # <--- Needed by variable substitution dict lookups

if __name__ == '__main__':
    test_AnyOrderODESystem___1()
    test_ODESystem___2()
