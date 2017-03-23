#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from scipy.integrate import quad
from scipy.optimize import brentq
from argparse import ArgumentParser, RawTextHelpFormatter

parser = ArgumentParser(description=
"""Plots the momentum occupation fraction for electrons as a function
of momentum for a single given number density and multiple
temperatures.  Inspired by Fig. 3.2 of Onno Pols' lecture notes on
stellar evolution and computed using his equations (3.28)--(3.30).

https://www.astro.ru.nl/~onnop/education/stev_utrecht_notes/""")

parser.add_argument('-T', '--temperature', type=float, nargs='+',
                    default=[2e6, 2e7],
                    help="list of temperatures, in kelvin")
parser.add_argument('-n', '--density', type=float, default=6e27,
                    help="number density, in units of cm^{-3}")
args = parser.parse_args()

h = 6.626069e-27
k = 1.38065e-16
m_e = 9.109382e-28
clight = 2.997902458e10

# eq (3.28)
def n_max(p):
    return 8.*np.pi/h**3*p**2

# eq (3.29)
def n_e(p, psi, T):
    x = 2.*m_e*k*T
    return 2./h**3/(np.exp(p**2/x-psi)+1.)*4.*np.pi*p**2

# eq (3.30)
def n_MB(p, T):
    x = 2.*m_e*k*T
    return np.exp(-p**2/x)/(np.pi*x)**1.5*4.*np.pi*p**2

def get_psi(n, T):
    pmax = 2.*m_e*clight
    return brentq(lambda z: quad(n_e, 0., pmax, args=(z, T))[0]-n, 0., 1e10)
           
pmax = 2.*m_e*clight
prange = np.linspace(0., 1e-17, 201)
n0 = args.density

pl.plot(prange/m_e/clight, n_max(prange)*m_e*clight/n0, 'k-.')
for T in args.temperature:
    psi = get_psi(n0, T)

    c, = pl.plot(prange/m_e/clight, n_e(prange, psi, T)*m_e*clight/n0)
    pl.plot(prange/m_e/clight, n0*n_MB(prange, T)*m_e*clight/n0, '--', color=c.get_color())

pl.xlabel('momentum ($m_ec$)')
pl.ylabel('fractional number density per unit momentum')
pl.show()
