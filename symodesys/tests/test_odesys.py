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

def test_SimpleFirstOrderODESystem___1():
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

def _mk_brith_death_system(n):
    b, d = map(sympy.symbols, ['b:'+str(n),'d:'str(n)])
    b, d = b[1:], d[:-1] # First isn't born, last doesn't die
    t = sympy.Symbol('t')
    x = [sympy.Function('x'+i)(t) for i in range(n)]
    births = [0]+[b[i]*x[i-1] for i in range(1:n)] # birth processes
    deaths = [d[i]*x[i] for i in range(n-1)]+[0] # death processes
    class BirthDeath(FirstOrderODESystem):
        indepv = t
        param_symbs = b+d
        f = OrderedDict([(x[i].diff(t), births[i]-deaths[i]) for i in range(n)])

    return BirthDeath()

def test__ODESysBase___1():
    pass
    # test mk_func


def test_FirstOrderODESystem___1():
    n = 3
    bd = _mk_brith_death_system(n)

    b, d = map(sympy.symbols, ['b:'+str(n),'d:'str(n)])
    b, d = b[1:], d[:-1] # First isn't born, last doesn't die
    t = sympy.Symbol('t')
    x = [sympy.Function('x'+i)(t) for i in range(n)]
    births = [0]+[b[i]*x[i-1] for i in range(1:n)] # birth processes
    deaths = [d[i]*x[i] for i in range(n-1)]+[0] # death processes

    exprs = [births[i]-deaths[i] for i in range(n)]
    # Test __getitem__
    assert bd['x0'] == x[0]
    assert bd['x1'] == x[1]
    assert bd['x2'] == x[2]
    assert bd.all_depv == x
    assert bd.non_analytic_depv == x
    assert bd.analytic_depv == []
    assert bd.solved_exprs == []
    assert bd.param_and_sol_symbs == b+d
    assert bd.known_symbs == [t]+x+b+d
    assert bd == _mk_brith_death_system(n)

    x_ = [sympy.Symbol('x'+i) for i in range(n)]
    assert bd.forbidden_symbs == [t]+x+b+d+x_

    y = [sympy.Function('y'+i)(t) for i in range(n)]
    bd.subs(dict(zip(x,y)))
    assert bd.is_first_order
    assert bd.get_highest_order() == 1
    assert bd.is_autonomous
    assert map(bd.unfunc_depv, x) == x_
    bd._do_sanity_check_of_odeqs()
    assert bd.eqs == [sympy.Eq(x[i].diff(t), expr[i]) for i in range(i)]


if __name__ == '__main__':
    test_AnyOrderODESystem___1()
    test_FirstOrderODESystem___1()
    test_SimpleFirstOrderODESystem___1()
