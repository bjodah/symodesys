#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
import os
import sys
import pickle
import logging
from functools import wraps
from collections import OrderedDict # for OrderedDefaultdict
from hashlib import md5


# other imports
import sympy
from sympy.utilities.autowrap import autowrap, ufuncify


def deprecated(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        print('This is a deprecation warning regarding: {}'.format(
            f.__name__))
        return f(*args, **kwargs)
    return wrapper


def cache(f):
    data = {}
    @wraps(f)
    def wrapper(*args):
        hashable_args=[]
        for x in args:
            if isinstance(x, dict):
                hashable_args.append(frozenset(x.items()))
            elif isinstance(x, list):
                hashable_args.append(tuple(x))
            else:
                hashable_args.append(x)
        hashable_args = tuple(hashable_args)
        if not hashable_args in data:
            data[hashable_args] = f(*args)
        return data[hashable_args]
    wrapper.cache_clear = lambda: data.clear()
    wrapper.cache = data
    return wrapper


def subs_set(s, subsd):
    """
    Substititues entities in a set ``s'' for matching keys in ``subsd''
    with corresponding values
    """
    t = s.copy()
    for k, v in subsd.iteritems():
        if k in s:
            t.remove(k)
            t.add(v)
    return t


def log_call_debug(func):
    from functools import wraps
    @wraps(func)
    def wrapper(*args, **kwargs):
        logging.debug(func.__name__ + ' called.')
        return func(*args, **kwargs)
    return wrapper


class OrderedDefaultdict(OrderedDict):
    """
    From http://stackoverflow.com/questions/4126348/\
    how-do-i-rewrite-this-function-to-implement-ordereddict/4127426#4127426
    """

    def __init__(self, *args, **kwargs):
        newdefault = None
        newargs = ()
        if args:
            newdefault = args[0]
            if not (newdefault is None or callable(newdefault)):
                raise TypeError('first argument must be callable or None')
            newargs = args[1:]
        self.default_factory = newdefault
        super(self.__class__, self).__init__(*newargs, **kwargs)

    def __missing__ (self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):  # optional, for pickle support
        args = self.default_factory if self.default_factory else tuple()
        return type(self), args, None, None, self.items()
