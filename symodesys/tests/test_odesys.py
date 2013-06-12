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
    l = sympy.symbols('l:'+str(n-1))+[0] # Death prob. coefficient
    x = [sympy.Function('x'+i)(t) for i in range(n)]
    d = [x[i]*l[i] for i in range(n)] # death rate
    b = [0]+[x[i-1]*l[i] for i in range(1:n)] # birth rate
    t = sympy.Symbol('t')
    class BirthDeath(FirstOrderODESystem):
        indepv = t
        param_symbs = b+d
        f = OrderedDict([(x[i].diff(t), b[i]-d[i]) for i in range(n)])

    return BirthDeath()

def test__ODESysBase___1():
    pass
    # test mk_func


def test_FirstOrderODESystem___1():
    n = 3
    bd = _mk_brith_death_system(n)

    l = sympy.symbols('l:'+str(n-1))+[0] # Death prob. coefficient
    x = [sympy.Function('x'+i)(t) for i in range(n)]
    d = [x[i]*l[i] for i in range(n)] # death rate
    b = [0]+[x[i-1]*l[i] for i in range(1:n)] # birth rate
    t = sympy.Symbol('t')
    exprs = [b[i]-d[i] for i in range(n)]
    # _ODESystemBase
    #   __getitem__
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

    # FirstOrderODESystem specific
    assert bd._odeqs == OrderedDict([
        (y[i], (1, exprs[i])) for i in range(n)])
    assert all([bd.f[k] == exprs[i] for i,k in enumerate(bd.all_depv)])

    # _ODESystemBase
    y = [sympy.Function('y'+i)(t) for i in range(n)]
    y_ = [sympy.Symbol('y'+i) for i in range(n)]
    bd.subs(dict(zip(x,y)))
    assert bd.is_first_order
    assert bd.get_highest_order() == 1
    assert bd.is_autonomous
    assert map(bd.unfunc_depv, y) == y_
    bd._do_sanity_check_of_odeqs()
    assert bd.eqs == [sympy.Eq(y[i].diff(t), expr[i]) for i in range(i)]
    assert bf.eq(bd['y0']) == sympy.Eq(y[0].diff(t), -d[0]*y[0])
    assert bd.from_list_of_eqs(bd.eqs) == bd
    C0 = sympy.Symbol('C0')
    assert bd.attempt_analytic_sol(y[0], C0*sympy.exp(-d[0]*t), [C0])
    C1 = sympy.Symbol('C1')
    assert not bd.attempt_analytic_sol(y[1], C1*sympy.exp(-d[1]*t), [C1])
    assert bd.is_linear

    bd.recursive_analytic_auto_sol()


if __name__ == '__main__':
    test_AnyOrderODESystem___1()
    test_FirstOrderODESystem___1()
    test_SimpleFirstOrderODESystem___1()
