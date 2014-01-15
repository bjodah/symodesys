#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import argparse

import numpy as np
import matplotlib.pyplot as plt
try:
    from cInterpol import PiecewisePolynomial
except ImportError:
    from scipy.interpolate import PiecewisePolynomial

def main(input, output, nplotx, nderiv):
    """
    """
    print(input)
    data = np.genfromtxt(input[0])
    ny = (data.shape[1]-1) // nderiv
    plotx = np.linspace(data[0,0], data[-1,0], nplotx)
    ls = ['-','--',':']
    color = ['r','g','b','k']
    for i in range(ny):
        pp = PiecewisePolynomial(np.ascontiguousarray(data[:,0]),
                                 np.ascontiguousarray(
            data[:,1+i*nderiv:1+(i+1)*nderiv]))
        plt.plot(plotx, pp(plotx), label=str(i),
                 ls=ls[i%len(ls)], color=color[i%len(color)])
    plt.legend()
    plt.savefig(output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type=str, default='plot.png')
    parser.add_argument('--nplotx', type=int, default=1000)
    parser.add_argument('--nderiv', type=int, default=3)
    parser.add_argument('input', nargs=1, type=str)
    argd = vars(parser.parse_args())
    main(**argd)
