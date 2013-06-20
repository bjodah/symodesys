#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
import os
import sys
import pickle
from functools import wraps
from collections import OrderedDict # for OrderedDefaultdict
from hashlib import md5


# other imports
import sympy
from sympy.utilities.autowrap import autowrap, ufuncify
from mako.template import Template
from mako.exceptions import text_error_template


def render_mako_template_to(template, outpath, subsd):
    """
    template: either string of path or file like obj.
    """
    if hasattr(template, 'read'):
        # set in-file handle to provided template
        ifh = template
    else:
        # Assume template is a string of the path to the template
        ifh = open(template, 'rt')

    template_str = ifh.read()
    with open(outpath, 'wt') as ofh:
        try:
            rendered = Template(template_str).render(**subsd)
        except:
            print(text_error_template().render())
            raise

        ofh.write(rendered)


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
    """
    Imports (cython generated) shared object file (.so)

    Warning, Python's lazy caching is horrible for reimporting
    same path of an .so file. It will not detect the new time stamp
    nor new checksum but will use old module. Use unique names for
    this reason
    """
    import imp
    path, name = os.path.split(filename)
    name, ext = os.path.splitext(name)
    fobj, filename, data = imp.find_module(name, [path])
    mod = imp.load_module(name, fobj, filename, data)
    return mod


def import___not_working(filename):
    """
    Imports shared object file (.so)
    """
    path, name = os.path.split(filename)
    name, ext = os.path.splitext(name)
    if not path in sys.path:
        sys.path.append(path)
    for k,v in sys.modules.items():
        if hasattr(v, '__file__'):
            if v.__file__.endswith('wrapper.so'):
                print k,v
    if name in sys.modules:
        #from IPython.extensions.autoreload import superreload
        #print 'superreloading: '+name
        #mod = superreload(sys.modules[name])
        # print 'pop in modules'
        # sys.modules.pop(name)
        if name in sys.meta_path:
            print 'pop in meta path'
            sys.meta_path.pop(name)
        #mod = __import__(name)
        #reload(mod)
        import importlib
        mod = importlib.import_module(name)
    else:
        mod = __import__(name)
    print 'id of .transform function: ', id(mod.transform)
    print 'md5 of .so file: ',md5_of_file(mod.__file__).hexdigest()
    return mod


def missing_or_other_newer(path, other_path):
    if not os.path.exists(path) or \
       os.path.getmtime(other_path) > os.path.getmtime(path):
        return True
    return False



class HasMetaData(object):
    """
    Provides convenice methods for a class to pickle some metadata
    """
    metadata_filename = '.metadata'

    @classmethod
    def _get_metadata_key(cls, kw):
        """ kw could be e.g. 'compiler' """
        return cls.__name__+'_'+kw


    @classmethod
    def get_from_metadata_file(cls, dirpath, key):
        """
        Get value of key in metadata file dict.
        """
        fullpath = os.path.join(dirpath, cls.metadata_filename)
        if os.path.exists(fullpath):
            d = pickle.load(open(fullpath,'r'))
            return d[key] #.get(cls._get_metadata_key(key), None)
        else:
            raise IOError("No such file: {}".format(fullpath))

    @classmethod
    def save_to_metadata_file(cls, dirpath, key, value):
        """
        Store `key: value` in metadata file dict.
        """
        fullpath = os.path.join(dirpath, cls.metadata_filename)
        if os.path.exists(fullpath):
            d = pickle.load(open(fullpath,'r'))
            d.update({key: value})
            pickle.dump(d, open(fullpath,'w'))
        else:
            pickle.dump({key: value}, open(fullpath,'w'))
