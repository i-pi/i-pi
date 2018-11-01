"""Decorators for easy utilizing i-PI simulations and more."""

# This file is part of i-PI.
# i-PI Copyright (C) 2014-2015 i-PI developers
# See the "licenses" directory for full license information.


def cached(f):
    """Cache decorator."""

    _cache = {}

    def func(*args):
        if args in _cache:
            return _cache[args]
        res = f(*args)
        _cache[args] = res
        return res

    return func
