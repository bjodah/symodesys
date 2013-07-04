#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division

import argparse

import numpy as np
import matplotlib.pyplot as plt
from cInterpol import PiecewisePolynomial

def main(input, output, nx=1000):
    """
    """
    print(input)
    data = np.genfromtxt(input[0])
    ny = (data.shape[1]-1) // 3
    plotx = np.linspace(data[0,0], data[-1,0], nx)
    ls = ['-','--',':']
    color = ['r','g','b','k']
    for i in range(ny):
        pp = PiecewisePolynomial(data[:,0], data[:,1+i*3:1+(i+1)*3])
        plt.plot(plotx, pp(plotx), label=str(i),
                 ls=ls[i%len(ls)], color=color[i%len(color)])
    plt.legend()
    plt.savefig(output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type=str, default='plot.png')
    parser.add_argument('-nx', type=int, default=1000)
    parser.add_argument('input', nargs=1, type=str)
    argd = vars(parser.parse_args())
    main(**argd)
