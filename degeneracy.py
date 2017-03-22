#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from scipy.integrate import quad

# from Onno Pols' notes, Sec. 3.3.5

h = 6.626069e-27
k = 1.38065e-16
m_e = 9.109382e-28
clight = 2.997902458e10

def n_e(p, psi, T):
    x = 2.*m_e*k*T
    return 2./h**3/(np.exp(p**2/x-psi)+1.)*4.*np.pi*p**2

def n_MB(p, T):
    x = 2.*m_e*k*T
    return np.exp(-p**2/x)/(np.pi*x)**1.5*4.*np.pi*p**2

def n_max(p):
    return 8.*np.pi/h**3*p**2

pmax = 2.*m_e*clight
prange = np.linspace(0., 1e-17, 201)
# for ne=6e27 at 2e5/6/7K: 698.6083467, 69.84917591, 6.8640177
psi = 69.84917591
T = 2e6
n = quad(n_e, 0., pmax, args=(psi, T))
print(n)
pl.plot(prange, n_e(prange, psi, T))
pl.plot(prange, n[0]*n_MB(prange, T))
pl.plot(prange, n_max(prange), 'k--')
pl.show()
