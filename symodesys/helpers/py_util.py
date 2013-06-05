#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
import os
import imp
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


def md5_of_file(path):
    """
    Use .digest() or .hexdigest() on returned object
    to get binary or hex encoded string.
    """
    md = md5()
    with open(path,'rb') as f:
        for chunk in iter(lambda: f.read(128*md.block_size), b''):
             md.update(chunk)
    return md


def subs_set(s, subsd):
    """
    Substititues entities in a set ``s'' for matching keys in ``subsd''
    with corresponding values
    """
    t = s.copy()
    for k, v in subsd.iteritems():
        if k in s:
            s.remove(k)
            s.add(v)
    return t


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




def import_(filename):
    """ Imports (cython generated) shared object file (.so) """
    path, name = os.path.split(filename)
    name, ext = os.path.splitext(name)
    fobj, filename, data = imp.find_module(name, [path])
    mod = imp.load_module(name, fobj, filename, data)
    return mod
