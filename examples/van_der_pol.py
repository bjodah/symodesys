#!/usr/bin/env python
# -*- coding: utf-8 -*-

from symodesys.firstorder import FirstOrderODESystem

class VanDerPolOscillator(FirstOrderODESystem):
    ny = 2
    nk = 1

    @property
    def f(self):
        u, v = *self.y
        mu = self.k[0]
        return [v,
                -u + mu * v * (1 - u ** 2),
                ]

