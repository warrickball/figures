#!/usr/bin/env python

"""Oscillations of an isothermal sphere, following formulae given by
Gough (2003; Ap&SS 284, 165), eq. (21):

https://link.springer.com/content/pdf/10.1023%2FA%3A1023275232004.pdf

L = \ell + 1/2
\Psi = \Psi_0 r^(1/2) J_L(j_{Ln}r/R)P_l^m(\cos\theta)\cos(m\phi)\cos(\omega t)

where

j_{Ln} is nth zero of J_L

"""

import numpy as np
from matplotlib import pyplot as pl
from scipy.optimize import brentq

# J(k,x) is Bessel function of first kind, kth order
from scipy.special import jv as J

# lpmv(m,l,mu) is associated Legendre poly of azimuthal order m and
# angular degree l
from scipy.special import lpmv

def get_J_zeros(l, z0s):
    zs = []
    for i in range(1, len(z0s)):
        try:
            zs.append(brentq(lambda z: J(0.5+l, z), z0s[i-1], z0s[i]))
        except ValueError:
            continue
            
    return np.array(zs)

Dz = np.pi
x = np.linspace(0., 24.*np.pi, 101)
z0 = np.arange(0.01, x[-1], np.pi/2)
for l in range(4):
    z = get_J_zeros(l, z0)
    pl.plot(np.mod(z/Dz-0.01, 1.0)+0.01, z/Dz, 'o', label=r'$\ell=%i$' % l)

pl.xlim(-0.05, 1.05)
pl.xlabel('frequency modulo large separation (=$\pi$)')
pl.ylabel('frequency/large separation ($\\approx$ radial order)')
pl.title('echelle diagram of an isothermal sphere')
pl.legend()
pl.show()
