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


def _make_hashable(x):
    if isinstance(x, dict):
        return frozenset(x.items())
    elif isinstance(x, list):
        return tuple(x)
    else:
        return x


def hashable(it):
    frozen = []
    for x in it:
        candidate = _make_hashable(x)

        try:
            iter_cand = iter(candidate)
            candidate = hashable(iter_cand)
        except TypeError:
            pass
        frozen.append(candidate)

    return tuple(frozen)


def cache(f):
    """
    Defines an infinte cache decorator for a function
    not accepting keyword arguments.
    """
    data = {}

    @wraps(f)
    def wrapper(*args):
        hashable_args=[]
        map(hashable_args.append, map(_make_hashable, args))
        hashable_args = tuple(hashable_args)
        if not hashable_args in data:
            data[hashable_args] = f(*args)
        return data[hashable_args]

    wrapper.cache_clear = lambda: data.clear()
    wrapper.cache = data
    return wrapper


def cache_w_kwargs(*kwtup):
    """
    Decorator factory.
    >>> @cache_w_kwargs('x','y')
    ... def func(x, y=None):
    ...     pass
    """
    def cch(f):
        data = {} # Could be made to size-limited LRU

        @wraps(f)
        def wrapper(*args, **kwargs):
            hshbl = [] # hashable
            map(hshbl.append, map(_make_hashable, args))
            for kw in kwtup[len(args):]:
                hshbl.append(_make_hashable(kwargs[kw]))
            hshbl = tuple(hshbl)
            if not hshbl in data:
                data[hshbl] = f(*args, **kwargs)
            return data[hshbl]
        wrapper.cache_clear = lambda: data.clear()
        wrapper.cache = data
        return wrapper
    return cch


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
