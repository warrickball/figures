#!/usr/bin/env python

"""Creates 2D plots of power in an oscillation mode multiplet, as in
Gizon & Solanki (2003)."""

import numpy as np
from matplotlib import pyplot as pl

# hardcoded up to l=4 for speed
energy_list = [
    [lambda x: 1],
    [lambda x: np.cos(x)**2, 
     lambda x: np.sin(x)**2/2.],
    [lambda x: (3.*np.cos(x)**2-1)**2/4.,
     lambda x: np.sin(2.*x)**2*3./4.,
     lambda x: np.sin(x)**4*3./4.],
    [lambda x: (5.*np.cos(3.*x)+3.*np.cos(x))**2/64.,
     lambda x: (5.*np.cos(2.*x)+3.)**2*3./64.,
     lambda x: np.cos(x)**2*np.sin(x)**4*15./8.,
     lambda x: np.sin(x)**6*5./16.],
    [lambda x: (35.*np.cos(x)**4-30.*np.cos(x)**2+3.)**2/64.,
     lambda x: (7./2.*np.sin(4.*x)+np.sin(2.*x))**2*5./256.,
     lambda x: (7.*np.cos(2.*x)+5.)**2*np.sin(x)**4*5./128.,
     lambda x: np.cos(x)**2*np.sin(x)**6*35./16.,
     lambda x: np.sin(x)**8*35./128.]]


def lorentz(x): return 1./(1.+x**2)

def multiplet(x, l, freq, width, splitting, inc):
    y = energy_list[l][0](inc)*lorentz((x-freq)/width)
    for m in range(1,l+1):
        y = y + energy_list[l][abs(m)](inc)*lorentz((x-freq-m*splitting)/width)
        y = y + energy_list[l][abs(m)](inc)*lorentz((x-freq+m*splitting)/width)

    return y

l = 2
x = np.linspace(-l-0.5,l+0.5,1000)    # in units of rotational splitting
inc = np.linspace(0., np.pi/2., 300)
y = np.vstack([multiplet(x, 2, 0., 0.2, 1., inci) for inci in inc])

pl.imshow(-y, origin='lower', cmap='gray', extent=[-l-0.5,l+0.5,0,90],
          aspect='auto')
pl.grid('off')
# pl.xlabel(r'$(\nu-\nu_\mathrm{mode})/\delta\nu_\mathrm{splitting}$')
pl.xlabel('frequency shift / rotation frequency')
pl.ylabel(r'inclination ($^\circ$)')
pl.show()
