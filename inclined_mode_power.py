#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from scipy.special import lpmn
from argparse import ArgumentParser

parser = ArgumentParser(description="Creates 2D plots of power in an oscillation "
                        "mode multiplet, as in Gizon & Solanki (2003).")
parser.add_argument('-l', '--ell', type=int, default=2,
                    help="angular degree (default=2)")
parser.add_argument('-w', '--width', type=float, default=0.2,
                    help="linewidth (default=0.2)")
args = parser.parse_args()

def factorial(n):
    if n == 0:
        return 1
    elif n > 0:
        return np.prod(range(1, n+1))
    else:
        raise ValueError('can only return factorial of integer >= 0')

def lorentz(x): return 1./(1.+x**2)

def multiplet(x, l, freq, width, splitting, inc):
    legendre2 = lpmn(l, l, np.cos(inc))[0]**2
    y = lorentz((x-freq)/width)*legendre2[0][l]
    for m in range(1, l+1):
        energy = legendre2[m][l]*factorial(l-m)/factorial(l+m)
        y = y + lorentz((x-freq-m*splitting)/width)*energy
        y = y + lorentz((x-freq+m*splitting)/width)*energy

    return y


l = args.ell
w = args.width

x = np.linspace(-l-0.5,l+0.5,1000)    # in units of rotational splitting
inc = np.linspace(0., np.pi/2., 300)
y = np.vstack([multiplet(x, l, 0., w, 1., inci) for inci in inc])

pl.imshow(-y, origin='lower', cmap='gray',
          extent=[-l-0.5, l+0.5, 0, 90], aspect='auto')
pl.grid(False)
# pl.xlabel(r'$(\nu-\nu_\mathrm{mode})/\delta\nu_\mathrm{splitting}$')
pl.xlabel('frequency shift / rotation frequency')
pl.ylabel(r'inclination ($^\circ$)')
pl.show()
