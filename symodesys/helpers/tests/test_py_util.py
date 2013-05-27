#!/usr/bin/env python
# -*- coding: utf-8 -*-

from ..py_util import subs_set

def test_subs_set():
    s1 = set('abc')
    d = dict('b':'d')
    s2 = subs_set(s, d)
    assert s2 == set('acd')

if __name__ == '__main__':
    test_subs_set()
