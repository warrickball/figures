#!/usr/bin/env python

"""Plots a representation of the radial displacement of a solar
oscillation mode using Mayavi (as opposed to matplotlib).

"""

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-l', '--ell', type=int, default=20,
                    help="Angular degree l (default l=20)")
parser.add_argument('-m', '--emm', type=int, default=16,
                    help="Azimuthal order m (default m=16)")
parser.add_argument('-n', '--enn', type=int, default=None,
                    help="Radial order n.  You can only choose one of this "
                    "or the cyclic frequency n (-f, --freq).  (default n=14)")
parser.add_argument('-f', '--freq', type=float, default=None,
                    help="Cyclic frequency in mHz.  You can only choose one "
                    "of this or the radial order n (-n, --enn).")
parser.add_argument('-o', '--output', type=str, default=None,
                    help="save figure to given filename without displaying "
                    "it (forces software rendering)")
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
from tomso import adipls, fgong

if args.freq and args.enn:
        raise ValueError("""you can only use one of the radial order (-n,
        --enn) or the cyclic frequency (-f, --freq)""")

if args.output:
        mlab.options.offscreen = True

mlab.figure(1, bgcolor=tuple(args.bgcolor), fgcolor=(0, 0, 0), size=(600, 600))
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

    S = fgong.load_fgong('data/modelS.fgong', return_object=True)
    amde = adipls.load_amde('data/modelS.amde', return_object=True)
    I = np.where(amde.l==l)[0]
    if args.freq:
        i = I[np.argmin((amde.nu_Ri[I]-args.freq/1e3)**2)]
    elif args.enn:
        i = I[amde.n[I]==args.enn][0]
    else:
        i = I[amde.n[I]==14][0]
        
    print('   n = %i' % amde.n[i])
    print('freq = %.6f mHz' % (amde.nu_Ri[i]*1e3))
    r = amde.eigs[i][:,0]
    y1 = amde.eigs[i][:,1]
    rho = np.interp(r, S.x[::-1], S.rho[::-1])

    # r = np.linspace(0.,1.,51)
    # s = np.outer(np.sin(1.5*np.pi*r), np.ones(len(th)))
    s = np.outer(r*rho**0.5*y1, np.ones(len(th)))
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

