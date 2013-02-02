from functools import wraps

def cache(f):
    data = {}
    @wraps(f)
    def wrapper(*args):
        if args in data:
            return data[args]
        else:
            data[args] = f(*args)
            return data[args]
    return wrapper
