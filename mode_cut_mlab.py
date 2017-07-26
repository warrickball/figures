#!/usr/bin/env python

"""Plots a representation of the radial displacement of a solar
oscillation mode using Mayavi (as opposed to matplotlib).

"""

import numpy as np
from matplotlib import pyplot as pl
from mayavi import mlab
from scipy.special import sph_harm, lpmv, factorial
from scipy.optimize import fsolve
from tomso import adipls, io

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(600, 600))
mlab.clf()
ell = 20
emm = 16

# ph = np.linspace(0.5*np.pi, 2.0*np.pi, 100)
@mlab.show
def myplot():
    th = np.linspace(0.0, np.pi, 201)
    ph = np.linspace(0.5*np.pi, 2.0*np.pi, 301)
    Th, Ph = np.meshgrid(th, ph)

    x = np.outer(np.cos(ph), np.sin(th))
    y = np.outer(np.sin(ph), np.sin(th))
    z = np.outer(np.ones(np.size(ph)), np.cos(th))
    s = sph_harm(emm, ell, Ph, Th).real

    mlab.mesh(x, y, z, scalars=s, colormap='seismic')

    glob, var = io.load_fgong('data/modelS.fgong')
    css, eigs = adipls.load_amde('data/modelS.amde')
    I = np.where(css['ell']==ell)[0]
    i = I[np.argmin((css['nu_Ri'][I]-3.)**2)]
    r = eigs[i][:,0]
    y1 = eigs[i][:,1]
    rho = np.interp(r, var[::-1,0]/glob[1],
                    var[::-1,4])

    # r = np.linspace(0.,1.,51)
    # s = np.outer(np.sin(1.5*np.pi*r), np.ones(len(th)))
    s = np.outer(r*rho**0.5*y1, np.ones(len(th)))
    s = s*sph_harm(emm, ell, 
                   np.outer(np.ones(len(r)), np.ones(len(th))), 
                   np.outer(np.ones(len(r)), th)).real
    # s = 0.5+0.5*s/np.max(np.abs(s))
    smax = np.max(np.abs(s))/2.  # divide by two gives richer colour
    smin = -smax  # guarantees symmetry in colourmap

    x = np.outer(r, np.sin(th))
    y = 0.*x
    z = np.outer(r, np.cos(th))
    mlab.mesh(x, y, z, scalars=s, colormap='seismic',
              vmin=smin, vmax=smax)

    y = np.outer(r, np.sin(th))
    x = 0.*y
    z = np.outer(r, np.cos(th))
    mlab.mesh(x, y, z, scalars=s, colormap='seismic',
              vmin=smin, vmax=smax)

    mlab.view(azimuth=15.0, elevation=90.0, distance=4.2)

myplot()
