#!/usr/bin/env python3
#
# Exact solution of 1D transverse-field Ising model with periodic boundary conditions

from math import sqrt, cos, pi

L = 32
Gamma = 1

energy = 0
for i in range(L):
    k = pi * (2 * i + 1) / L
    energy += -sqrt(1 + Gamma**2 + 2 * Gamma * cos(k))

print(energy)
