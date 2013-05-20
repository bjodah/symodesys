#!/usr/bin/env python
# -*- coding: utf-8 -*-

# stdlib imports
from math import floor, sqrt, log
from itertools import chain


isqrt = lambda x: int(floor(sqrt(x)))
ilog2 = lambda x: int(log(x)/log(2))

class PrimeCache(object):
    def __init__(self):
        self._m = 2
        self.primes = set([2])

    def _explore(self, x):
        for i in xrange(self._m, x+1):
            for p in chain(self.primes, xrange(max(self.primes)+2, isqrt(i), 2)):
                if i % p == 0:
                    # i not a prime
                    break
            else:
                # i is a prime
                self.primes.add(i)
        self._m = max(self._m, x)

    def is_prime(self, x):
        self._explore(x)
        return x in self.primes

    def factorize(self, n):
        # factorizes number n
        i = isqrt(n)
        self._explore(i) # ensure primes explored up to i
        fctrs = []
        for p in self.primes:
            if i % p == 0:
                return [p] + self.factorize(i/p)
            if p > i:
                break
        self.primes.add(n)
        return [n]

class SortedPrimes(list):
    def __init__(self, start):
        self.x = start
        self._pc = PrimeCache()

    def __iter__(self):
        while True:
            if self._pc.is_prime(self.x): yield self.x
            self.x += 1
        raise StopIteration


def test_SortedPrimes():
    primes = []
    for i in SortedPrimes(17):
        primes.append(i)
        if i >= 101: break
    assert primes == [17, 19, 23, 29, 31, 37, 41, 43, 47, 53,
                     59, 61, 67, 71, 73, 79, 83, 89, 97, 101]

def main(start, stop):
    """
    Example of using primes
    """
    for i in SortedPrimes(start):
        print i
        if i >= stop: break


if __name__ == '__main__':
    main(17, 101)
