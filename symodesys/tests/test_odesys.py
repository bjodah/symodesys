import sympy

from collections import OrderedDict

from symodesys.odesys import _ODESystemBase, FirstOrderODESystem, SimpleFirstOrderODESystem, AnyOrderODESystem



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

def _mk_decay_chain_system(n):
    l = list(sympy.symbols('l:'+str(n-1))) # Decay prob. coefficient
    t = sympy.Symbol('t')
    x = [sympy.Function('x'+str(i))(t) for i in range(n)]
    d = [x[i]*l[i] for i in range(n-1)]+[0] # death rate
    b = [0]+[x[i-1]*(l+[0])[i] for i in range(1,n)] # birth rate
    class DecayChain(FirstOrderODESystem):
        indepv = t
        param_symbs = l
        f = OrderedDict([(x[i], b[i]-d[i]) for i in range(n)])

    return DecayChain()

def test__ODESysBase___1():
    pass
    # test mk_func


def test_FirstOrderODESystem___1():
    n = 3
    dc = _mk_decay_chain_system(n)

    l = list(sympy.symbols('l:'+str(n-1))) # Decay prob. coefficient
    t = sympy.Symbol('t')
    x = [sympy.Function('x'+str(i))(t) for i in range(n)]
    d = [x[i]*l[i] for i in range(n-1)]+[0] # death rate
    b = [0]+[x[i-1]*(l+[0])[i] for i in range(1,n)] # birth rate
    exprs = [b[i]-d[i] for i in range(n)]

    # For 3 coupled decays we can solve by hand and compare
    # solution:

    I = sympy.symbols('I:3')
    e = sympy.exp

    # Use integrating factor to solve by hand:
    analytic = {x[0]: I[0]*e(-l[0]*t),
                x[1]: l[0]*I[0]*(e(-l[0]*t)-e(-l[1]*t))/(l[1]-l[0])+I[1]*e(-l[1]*t), # l[0] != l[1]
                x[2]: l[1]*I[0]/(l[1]-l[0])*(1-e(-l[0]*t))+\
                (l[0]*I[0]/(l[1]-l[0])-I[1])*e(-l[1]*t)+I[2]+I[1]
    }

    # _ODESystemBase
    #   __getitem__
    assert dc['x0'] == x[0]
    assert dc['x1'] == x[1]
    assert dc['x2'] == x[2]
    assert dc.all_depv == x
    assert dc.non_analytic_depv == x
    assert dc.analytic_depv == []
    assert dc.analytic_relations == []
    assert dc.param_symbs == l
    assert dc.analytic_sol_symbs == []
    assert dc.param_and_sol_symbs == l
    assert dc.known_symbs == [t]+x+l
    assert dc == _mk_decay_chain_system(n)
    assert dc.eqs == [sympy.Eq(x[i].diff(t), exprs[i]) for i in range(n)]

    x_ = [sympy.Symbol('x'+str(i), real=True) for i in range(n)]
    assert dc.forbidden_symbs == [t]+x+l+x_

    # FirstOrderODESystem specific
    assert dc._odeqs == OrderedDict([
        (x[i], (1, exprs[i])) for i in range(n)])
    assert all([dc.f[k] == exprs[i] for i,k in enumerate(dc.all_depv)])

    # _ODESystemBase
    y  = [sympy.Function('y'+str(i))(t) for i in range(n)]
    y_ = [sympy.Symbol('y'+str(i), real=True) for i in range(n)]
    dc = dc.transform_depv(dict(zip(y,x)), dict(zip(x,y)))
    assert dc.is_first_order
    assert dc.get_highest_order() == 1
    assert dc.is_autonomous
    assert map(dc.unfunc_depv, y) == y_
    dc._do_sanity_check_of_odeqs()
    assert dc.eq(dc['y0']) == sympy.Eq(y[0].diff(t), -d[0]*y[0])
    assert dc.from_list_of_eqs(dc.eqs) == dc
    C0 = sympy.Symbol('C0')
    assert dc.attempt_analytic_sol(y[0], C0*sympy.exp(-d[0]*t), [C0])
    C1 = sympy.Symbol('C1')
    assert not dc.attempt_analytic_sol(y[1], C1*sympy.exp(-d[1]*t), [C1])

    # FirstOrderODESystem specific
    assert dc.is_linear

def test_FirstOrderODESystem___2():
    """ Tests recursive_analytic_auto_sol """
    # recursive_analytic_auto_sol
    dc = _mk_decay_chain_system(3)
    dc.recursive_analytic_auto_sol()
    assert len(analytic_depv) == 3
    assert len(dc.non_analytic_depv) == 0
    assert len(dc.param_and_sol_symbs) == 6
    # let's see that our analytic solutions match
    for t in [1.0, 5.0, 17.0]:
        for l in ([3.0, 5.0, 7.0], [17.0, 11.0, 13.0]):
            for i, analytic_expr in enumerate(bd.analytic):
                assert sympy.N(bd.solved_exprs[i].subs(subsd_t)-\
                               bd.solved_exprs[i].subs(subsd_t0)) -\
                    sympy.N(analytic_expr.subs(subsd)) < 1e-10 # some tolerance

if __name__ == '__main__':
    test_AnyOrderODESystem___1()
    test_FirstOrderODESystem___1()
    test_SimpleFirstOrderODESystem___1()
    test_SimpleFirstOrderODESystem___2()
