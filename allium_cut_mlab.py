#!/usr/bin/env python

"""Plots a representation of the radial displacement of an oscillation
mode in a sphere with piecewise-constant sound speed using Mayavi (as
opposed to matplotlib).

"""

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-l', '--ell', type=int, default=20,
                    help="angular degree l (default l=20)")
parser.add_argument('-m', '--emm', type=int, default=16,
                    help="azimuthal order m (default m=16)")
parser.add_argument('-f', '--freq', type=float, default=1.0,
                    help="frequency in units of Δω (default=1.0)")
parser.add_argument('-c', type=float, nargs='+', default=[1.0],
                    help="sound speeds of each shell (default=[1.0])")
parser.add_argument('-R', type=float, nargs='+', default=[1.0],
                    help="outer radius of each shell (default=[1.0])")
parser.add_argument('-o', '--output', type=str, default=None,
                    help="save figure to given filename without displaying "
                    "it (forces software rendering)")
parser.add_argument('--resolution', type=float, nargs=2, default=[600,600],
                    help="resolution of image (default=[600,600])")
parser.add_argument('--view', type=float, nargs=2, default=[15.0, 90.0],
                    help="viewing angle (default=[15.0, 90.0])")
parser.add_argument('--bgcolor', type=float, nargs=3, default=[1,1,1],
                    help="background colour, as [0..1] RGB values "
                    "(default=1,1,1)")
parser.add_argument('--transparent', action='store_true',
                    help="save image to file with a transparent background")
args = parser.parse_args()

import pyface.qt
import numpy as np
from mayavi import mlab
from scipy.special import sph_harm, lpmv, factorial
from scipy.optimize import fsolve
from allium import Sphere

if args.output:
        mlab.options.offscreen = True

mlab.figure(1, bgcolor=tuple(args.bgcolor), fgcolor=(0, 0, 0), size=args.resolution)
mlab.clf()
l = args.ell
m = args.emm

# ph = np.linspace(0.5*np.pi, 2.0*np.pi, 100)
@mlab.show
def myplot():
    th = np.linspace(0.0, np.pi, 201)
    ph = np.linspace(0.5*np.pi, 2.0*np.pi, 301)
    Th, Ph = np.meshgrid(th, ph)

    x = np.outer(np.cos(ph), np.sin(th))
    y = np.outer(np.sin(ph), np.sin(th))
    z = np.outer(np.ones(np.size(ph)), np.cos(th))
    s = sph_harm(m, l, Ph, Th).real

    mlab.mesh(x, y, z, scalars=s, colormap='seismic')

    S = Sphere(c=args.c, R=args.R)
    S.search_eigenfrequencies(l, np.linspace(args.freq-2, args.freq+2, 11)*S.Delta_omega)

    i = np.argmin((np.array(S.eigenfrequencies[l])-args.freq*S.Delta_omega)**2)
    r = np.linspace(0, np.max(args.R), 1001)
    y = S.eigenfunction(l, i)(r)

    s = np.outer(y, np.ones(len(th)))
    s = s*sph_harm(m, l,
                   np.outer(np.ones(len(r)), np.ones(len(th))), 
                   np.outer(np.ones(len(r)), th)).real
    # s = 0.5+0.5*s/np.max(np.abs(s))
    smax = np.max(np.abs(s))/2.  # divide by two gives richer colour
    smin = -smax  # guarantees symmetry in colourmap

    # first inner semicircle
    x = np.outer(r, np.sin(th))
    y = 0.*x
    z = np.outer(r, np.cos(th))
    mlab.mesh(x, y, z, scalars=-s, colormap='seismic',
              vmin=smin, vmax=smax)

    # second inner semicircle
    y = np.outer(r, np.sin(th))
    x = 0.*y
    z = np.outer(r, np.cos(th))
    mlab.mesh(x, y, z, scalars=s, colormap='seismic',
              vmin=smin, vmax=smax)

    mlab.view(azimuth=args.view[0], elevation=args.view[1], distance=4.2)

myplot()
if args.output:
    if args.transparent:
        import matplotlib.pyplot as pl

        pl.imsave(args.output, mlab.screenshot(mode='rgba', antialiased=True))
    else:
        mlab.savefig(args.output)

