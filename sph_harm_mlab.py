#!/usr/bin/env python

import numpy as np
from mayavi import mlab
from scipy.special import sph_harm, lpmv
from scipy.optimize import fsolve
from numpy import pi, sin, cos
from argparse import ArgumentParser
import os

parser = ArgumentParser(description="""Uses Mayavi to plot a static spherical harmonic with a chosen
angular degree and azimuthal order.  """)
parser.add_argument('-l', '--ell', type=int, default=6,
                    help="angular degree")
parser.add_argument('-m', '--emm', type=int, default=3,
                    help="azimuthal order")
parser.add_argument('--save', type=str, default=None,
                    help="save plot to this file")
parser.add_argument('--Ntheta', type=int, default=101,
                    help="number of points in theta (latitude)")
parser.add_argument('--Nphi', type=int, default=101,
                    help="number of points in phi (longitude)")
parser.add_argument('-a', '--amplitude', type=float, default=1.0,
                    help="amplitude of oscillation")
parser.add_argument('--resolution', type=float, nargs=2, default=[400,400],
                    help="resolution of image")
parser.add_argument('--view', type=float, nargs=2, default=[45.0, 54.735610317245346],
                    help="viewing angle")
parser.add_argument('--distance', type=float, default=5.0,
                    help="camera distance")
parser.add_argument('--show-nodal-lines', dest='nodal_lines', action='store_true')
parser.add_argument('--hide-nodal-lines', dest='nodal_lines', action='store_false')
parser.set_defaults(nodal_lines=False)
args = parser.parse_args()

ell = args.ell
emm = args.emm

# Create a sphere
th = np.linspace(0., pi, args.Ntheta)
ph = np.linspace(-pi, pi, args.Nphi)
Th, Ph = np.meshgrid(th, ph)

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=args.resolution)
mlab.clf()
x = sin(Th)*cos(Ph)
y = sin(Th)*sin(Ph)
z = cos(Th)
s = sph_harm(emm,ell,Ph,Th).real
m = mlab.mesh(x, y, z, scalars=s, colormap='seismic')

# plot nodal lines
if args.nodal_lines:
# Get roots of assoc. Legendre polynomials

# this seems to work reasonably well. we basically just find zeros for
# too many initial guesses, then take the unique solutions.  some of
# these aren't zeros (because the root-finder fails) so we discard
# them.
    Nroots = ell - np.abs(emm) + 2
    mu = np.cos(np.linspace(0., np.pi, 5*Nroots))
    mu = np.squeeze([fsolve(lambda z: lpmv(emm, ell, z), mui) for mui in mu])
    mu = np.unique(np.around(mu, decimals=13))
    mu = np.array([mui for mui in mu
                   if np.around(lpmv(emm, ell, mui), decimals=4)==0.
                   and np.around(np.abs(mui), decimals=4) != 1.])

    node_kw = {'color':(0.,0.,0.),
               'line_width': 0.01, 'tube_radius': 0.01}
               # 'representation':'wireframe'}

    r = 1.001
    # equatorial
    for mui in mu:
        x = r*np.sqrt(1.-mui**2)*cos(ph)
        y = r*np.sqrt(1.-mui**2)*sin(ph)
        z = r*mui*np.ones(len(th))
           
        mlab.plot3d(x, y, z, **node_kw)

    # pole-to-pole
    for j in range(emm):
        Phi0 = pi*(2*j+1)/2/emm
        x = r*sin(th)*cos(Phi0)
        y = r*sin(th)*sin(Phi0)
        z = r*cos(th)

        mlab.plot3d(x, y, z, **node_kw)

        Phi0 = pi*(2*j+1)/2/emm + np.pi
        x = r*sin(th)*cos(Phi0)
        y = r*sin(th)*sin(Phi0)
        z = r*cos(th)

        mlab.plot3d(x, y, z, **node_kw)

# defaults are (45.0, 54.73561031724535, 6.744041908326433, array([0.0, 0.0, 0.0]))
mlab.view(azimuth=args.view[0], elevation=args.view[1], distance=args.distance)
if args.save:
    mlab.savefig(args.save)
    
mlab.show()
