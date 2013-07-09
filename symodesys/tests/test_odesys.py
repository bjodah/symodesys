import sympy
import sys
from collections import OrderedDict
from itertools import product

import numpy as np

from symodesys.odesys import (
    _ODESystemBase, FirstOrderODESystem, SimpleFirstOrderODESystem, AnyOrderODESystem)


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

    new_odesys = odesys.reduce_to_sys_of_first_order()

    #sympy.pprint(new_odesys.eqs)
    # TODO: add assert statements


def test_SimpleFirstOrderODESystem___1():
    class Decay(SimpleFirstOrderODESystem):
        real = True
        depv_tokens = 'u',
        param_tokens   = 'lambda_u',

        @property
        def expressions(self):
            return {self['u']: self['lambda_u'] * -self['u']}

    d=Decay()
    s = d.mk_depv('s')
    assert s.is_real
    assert d.mk_depv('s') == s # <--- Needed by variable substitution dict lookups


def _mk_decay_chain_system(n):
    l = list(sympy.symbols('l:'+str(n-1))) # Decay prob. coefficient
    t = sympy.Symbol('t')
    x = [sympy.Function('x'+str(i))(t) for i in range(n)]
    d = [x[i]*l[i] for i in range(n-1)]+[0] # death rate
    b = [0]+[x[i-1]*l[i-1] for i in range(1,n)] # birth rate
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
    b = [0]+[x[i-1]*l[i-1] for i in range(1,n)] # birth rate
    exprs = [b[i]-d[i] for i in range(n)]

    assert dc['x0'] == x[0]
    assert dc['x1'] == x[1]
    assert dc['x2'] == x[2]
    assert dc.all_depv == x
    assert dc.na_depv == x
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

    assert dc._odeqs == OrderedDict([
        (x[i], (1, exprs[i])) for i in range(n)])
    assert all([dc.f[k] == exprs[i] for i,k in enumerate(dc.all_depv)])

    # Numeric evaluation
    num_t0 = 13.0
    num_x0 = [2.0, 3.0, 5.0]
    num_l = [7.0, 11.0]

    assert np.allclose(dc.evaluate_f(num_t0, num_x0, num_l), np.array(
        [-2.0*7.0, 2.0*7.0 - 3.0*11.0, 3.0*11.0]))
    assert np.allclose(dc.evaluate_na_f(num_t0, num_x0, num_l), np.array(
        [-2.0*7.0, 2.0*7.0 - 3.0*11.0, 3.0*11.0]))

    ref_jac = np.array(
        [[-7.0,   0.0, 0.0],
         [ 7.0, -11.0, 0.0],
         [ 0.0,  11.0, 0.0]])
    assert np.allclose(dc.evaluate_jac(num_t0, num_x0, num_l), ref_jac)
    assert np.allclose(dc.evaluate_na_jac(num_t0, num_x0, num_l), ref_jac)

    assert np.allclose(dc.evaluate_dfdt(num_t0, num_x0, num_l), np.zeros(3))
    assert np.allclose(dc.evaluate_na_dfdt(num_t0, num_x0, num_l), np.zeros(3))

    ref_d2ydt2 = np.array(
        [-7.0*-2.0*7.0,
         7.0*-2.0*7.0 - 11.0*(2.0*7.0 - 3.0*11.0),
         11.0*(2.0*7.0 - 3.0*11.0)])
    assert np.allclose(dc.evaluate_d2ydt2(num_t0, num_x0, num_l), ref_d2ydt2)
    assert np.allclose(dc.evaluate_na_d2ydt2(num_t0, num_x0, num_l), ref_d2ydt2)

    y  = [sympy.Function('y'+str(i))(t) for i in range(n)]
    y_ = [sympy.Symbol('y'+str(i), real=True) for i in range(n)]
    dcy = dc.transform_depv(OrderedDict(zip(y,x)),
                            OrderedDict(zip(x,y)))
    dy = [y[i]*l[i] for i in range(n-1)]+[0] # death rate
    by = [0]+[y[i-1]*l[i-1] for i in range(1,n)] # birth rate

    assert dcy.f[dcy['y0']] == -l[0]*y[0]
    assert dcy.eq(dcy['y0']) == sympy.Eq(y[0].diff(t), -l[0]*y[0])
    assert dcy.is_first_order
    assert dcy.get_highest_order() == 1
    assert dcy.is_autonomous
    assert map(dcy.unfunc_depv, y) == y_
    dcy._do_sanity_check_of_odeqs()
    assert dcy.from_list_of_eqs(dcy.eqs).f == dcy.f

    assert dcy.is_linear
    ref_jac = sympy.Matrix(
        3,3, lambda i,j: (by[i]-dy[i]).diff(dcy.na_depv[j]))
    assert dcy.jac == ref_jac
    assert dcy.na_jac == ref_jac
    assert dcy.na_f == dcy.f
    assert dcy.na_dfdt == dcy.dfdt
    C0 = sympy.Symbol('C0')
    # attempt_analytic_sol is NOT idempotent
    assert dcy.attempt_analytic_sol(y[0], C0*sympy.exp(-l[0]*t), [C0])
    assert dcy.na_f != dcy.f
    assert dcy.na_jac != ref_jac
    assert dcy.na_dfdt != dcy.dfdt

    C1 = sympy.Symbol('C1')
    assert not dcy.attempt_analytic_sol(y[1], C1*sympy.exp(-l[1]*t), [C1])



def test_FirstOrderODESystem___2():
    """ Tests recursive_analytic_auto_sol """
    # recursive_analytic_auto_sol
    dc = _mk_decay_chain_system(3)
    dc.recursive_analytic_auto_sol()

    assert len(dc.analytic_depv) == 3
    assert len(dc.na_depv) == 0
    assert len(dc.param_and_sol_symbs) == 5

    l = list(sympy.symbols('l:2')) # Decay prob. coefficient

    undef_ref = [(l[1],l[0]), (l[1],0)] # what about l[0],0 ??
    undef = [sorted((eq.rhs,eq.lhs)) for eq in dc._solved_undefined]
    for a, b in undef_ref:
        assert sorted((a,b)) in undef

    t = sympy.Symbol('t')

    # For 3 coupled decays we can solve by hand and compare
    # solution:
    I = sympy.symbols('I:3')
    e = sympy.exp

    # Use integrating factor to solve by hand:
    analytic = [
        I[0]*e(-l[0]*t),
        l[0]*I[0]*(e(-l[0]*t)-e(-l[1]*t))/(l[1]-l[0])+\
            I[1]*e(-l[1]*t), # l[0] != l[1]
        I[0]/(l[1]-l[0])*(l[1]-l[0]+l[0]*e(-l[1]*t)-l[1]*e(-l[0]*t))+\
            I[2]+I[1]*(1-e(-l[1]*t)),
    ]

    # let's see that our analytic solutions match
    # symodesys solved expressoins for some (3x2x2)
    # combintaions of variable values of t, l and I
    ts = (1.0, 5.0, 17.0)
    ls = ([3.0, 5.0], [17.0, 11.0])
    Is = ([23.0, 29.0, 31.0], [17.0, 3.0, 43.0])
    for t_, l_, I_ in product(ts, ls, Is):
        subsd_t = {
            t: t_,
            l[0]: l_[0],
            l[1]: l_[1],
            I[0]: I_[0],
            I[1]: I_[1],
            I[2]: I_[2],
        }
        subsd_t0 = subsd_t.copy()
        subsd_t0.update({t: 0.0})

        # We need to solve sol symbs for t0:
        for depv, (expr, sol_symbs) in dc._solved.items():
            assert len(sol_symbs) == 1
            sol_symb = sol_symbs[0]
            j = dc.analytic_depv.index(depv)
            sols = sympy.solve(expr.subs(subsd_t0) - I_[j],
                                  sol_symb)
            assert len(sols) == 1
            sol_val = sols[0]
            subsd_t0.update({sol_symb: sol_val})
            subsd_t.update({sol_symb: sol_val})

        for i, analytic_rel in enumerate(analytic):
            # assert same within some tolerance
            assert sympy.N(dc.analytic_relations[i].subs(subsd_t)) -\
                sympy.N(analytic_rel.subs(subsd_t)) < 1e-10


if __name__ == '__main__':
    test_AnyOrderODESystem___1()
    test_FirstOrderODESystem___1()
    test_SimpleFirstOrderODESystem___1()

    skip_slow = False
    if len(sys.argv) > 1:
        if sys.argv[1] == '-s': # skip slow
            skip_slow = True

    if not skip_slow: test_FirstOrderODESystem___2()
