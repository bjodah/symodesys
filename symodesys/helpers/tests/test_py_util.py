# -*- coding: utf-8 -*-

from symodesys.helpers.py_util import subs_set, cache

def test_subs_set():
    s1 = set('abc')
    d = {'b': 'd'}
    s2 = subs_set(s1, d)
    assert s2 == set('acd')

def test_cache():
    global a
    a = 0

    @cache
    def f(x):
        global a
        a += 1
        return x**2

    assert a == 0
    assert f(2) == 4
    assert a == 1
    assert f(2) == 4
    assert a == 1
    assert f(3) == 9
    assert a == 2
