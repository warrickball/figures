#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as pl
from scipy.special import jn, jn_zeros
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-n', type=int, default=2,
                    help="radial order (default=2)")
parser.add_argument('-m', type=int, default=3,
                    help="angular/azimuthal order (default=3)")
parser.add_argument('--figsize', type=float, nargs=2,
                    help="figure size, passed to rcParams['figure.figsize'] "
                    "(you probably want to make it square)")
parser.add_argument('--levels', type=int, default=100,
                    help="number of levels passed to contourf (default 100)")
parser.add_argument('--padding', type=float, default=0.01,
                    help="fractional padding between edge and circle (default=0.01")
args = parser.parse_args()

if args.figsize:
    pl.rcParams['figure.figsize'] = args.figsize

# largely based on
# https://scipython.com/book/chapter-8-scipy/examples/drum-vibrations-with-bessel-functions/
# z(r,θ;t)=A Jm(kr) sin(mθ) cos(kνt)
# where k is nth zero of Jm
# by analogy with n,l,m of stellar oscillations

TAU = 2.*np.pi

# makes the formulae a bit neater
n = args.n # radial order (order of zero of bessel function)
m = args.m # angular order
k = jn_zeros(m, n+1)[n]

# r = np.linspace(0., 1., 51)
# th = np.linspace(0., 2.*np.pi, 51)
r, th = np.mgrid[0:1:101j, 0:1:101j]
th = th*TAU

x = r*np.cos(th)
y = r*np.sin(th)

k = jn_zeros(m, n+1)
if n > 0:
    z = jn(m, k[n]*r)*np.sin(m*th)
else:
    z = jn(m, k[n]*r)

scale = np.max(np.abs(z))
pl.clf()
pl.subplots_adjust(0, 0, 1, 1)
pl.contourf(x, y, z, levels=args.levels, cmap='seismic',
            vmin=-scale, vmax=scale)

th = np.linspace(0, TAU, 1001)
x = np.cos(th)
y = np.sin(th)
pl.plot(x, y, 'k-')

for ki in k[:-1]:
    r = ki/k[n]
    pl.plot(r*x, r*y, 'k--')

r = np.linspace(-1, 1, 101)
th = TAU/2/m if n > 0 else 0.
for i in range(m):
    pl.plot(r*np.cos(i*th), r*np.sin(i*th), 'k--')

b = 1.0 + args.padding
pl.axis([-b, b, -b, b])
pl.axis('off')

pl.show()
