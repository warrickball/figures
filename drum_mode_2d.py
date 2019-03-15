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
                    help="figure size, passed to rcParams['figure.figsize']")
parser.add_argument('--levels', type=int, default=100,
                    help="number of levels passed to contourf (default 100)")
parser.add_argument('--padding', type=float, default=0.01,
                    help="fractional padding between edge and circle (default=0.01")
parser.add_argument('-o', '--output', type=str,
                    help="save plot to file")
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

r, th = np.mgrid[0:1:201j, 0:1:201j]
th = th*TAU

x = r*np.cos(th)
y = r*np.sin(th)

k = jn_zeros(m, n+1)
if m > 0:
    z = jn(m, k[n]*r)*np.sin(m*th)
else:
    z = jn(m, k[n]*r)

zmax = np.max(np.abs(z))

pl.clf()
ax = pl.subplot(111, projection='polar')
b = args.padding
pl.subplots_adjust(top=1-b, bottom=b, left=b, right=1-b)
ax.contourf(th, r, z, levels=args.levels, cmap='seismic',
            vmin=-zmax, vmax=zmax)

th = np.linspace(0, TAU, 1001)
for ki in k[:-1]:
    r = ki/k[n]*np.ones_like(th)
    pl.plot(th, r, 'k--')

r = np.linspace(0, 1, 101)
th = TAU/2/m if m > 0 else 0.
for i in range(m):
    pl.plot(i*th*np.ones_like(r), r, 'k--')
    pl.plot(i*th*np.ones_like(r) + TAU/2, r, 'k--')

ax.set_xticklabels([])
ax.set_yticklabels([])

if args.output:
    pl.savefig(args.output)
else:
    pl.show()
