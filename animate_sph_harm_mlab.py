#!/usr/bin/env python

import numpy as np
from mayavi import mlab
from scipy.special import sph_harm, lpmv
from scipy.optimize import fsolve
from numpy import pi, sin, cos
from argparse import ArgumentParser
import os

parser = ArgumentParser(description="""Uses Mayavi to animates a spherical harmonic with a chosen
angular degree and azimuthal order.  """)
parser.add_argument('-l', '--ell', type=int, default=6,
                    help="angular degree")
parser.add_argument('-m', '--emm', type=int, default=3,
                    help="azimuthal order")
parser.add_argument('-o', '--output', type=str, default=None,
                    help="save figure to given filename without displaying "
                    "it (forces software rendering)")
parser.add_argument('--Ntheta', type=int, default=101,
                    help="number of points in theta (latitude)")
parser.add_argument('--Nphi', type=int, default=101,
                    help="number of points in phi (longitude)")
parser.add_argument('-a', '--amplitude', type=float, default=1.0,
                    help="amplitude of oscillation")
parser.add_argument('-P', '--period', type=float, default=1.0,
                    help="period of oscillation, in seconds")
parser.add_argument('--Nframes', type=int, default=40,
                    help="number of frames per oscillation")
parser.add_argument('--resolution', type=float, nargs=2, default=[400,400],
                    help="resolution of image")
parser.add_argument('--view', type=float, nargs=2, default=[45.0, 54.735610317245346],
                    help="viewing angle (default=45.0, 54.74)")
parser.add_argument('--bgcolor', type=float, nargs=3, default=[1,1,1],
                    help="background colour, as [0..1] RGB values "
                    "(default=1,1,1)")
parser.add_argument('--show-nodal-lines', dest='nodal_lines', action='store_true')
parser.add_argument('--hide-nodal-lines', dest='nodal_lines', action='store_false')
parser.set_defaults(nodal_lines=False)
args = parser.parse_args()

if args.output:
    mlab.options.offscreen = True
        
Nframes = args.Nframes
interval = int(args.period/Nframes*1000.)  # in milliseconds?
dphase = 2.*pi/Nframes
ell = args.ell
emm = args.emm

# Create a sphere
th = np.linspace(0., pi, args.Ntheta)
ph = np.linspace(-pi, pi, args.Nphi)
Th, Ph = np.meshgrid(th, ph)

x = sin(Th)*cos(Ph)
y = sin(Th)*sin(Ph)
z = cos(Th)
s = sph_harm(emm, ell, Ph, Th).real

mlab.figure(1, bgcolor=tuple(args.bgcolor), fgcolor=(0, 0, 0), size=args.resolution)
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
mlab.view(azimuth=args.view[0], elevation=args.view[1], distance=5.0)
# print(mlab.view())

# following http://zulko.github.io/blog/2014/11/29/data-animations-with-python-and-moviepy/
# duration = 1.0
def save_frame(filename, phase):
    dr = s*sin(phase)*args.amplitude
    x = (1.+dr)*sin(Th)*cos(Ph)
    y = (1.+dr)*sin(Th)*sin(Ph)
    z = (1.+dr)*cos(Th)

    m.mlab_source.set(x=x, y=y, z=z)
    # return mlab.screenshot(mode='rgb', antialiased=True)
    mlab.savefig(filename)

# # a = mpy.VideoClip(make_frame, duration=duration)
# # a.write_gif('sph_harm.gif', fps=20)
# dataset = [make_frame(t).T for t in np.arange(0.,1.,0.25)]
# print(len(dataset))
# print(dataset[0].shape)
# # mlab.show()
# write_gif(dataset, 'sph_harm.gif')
# save = True
if args.output:
    from subprocess import call
    phases = np.linspace(0., 2.*np.pi, Nframes+1)[:-1]
    for i, phase in enumerate(phases):
        save_frame('tmp/frame_%05i.png' % i, phase)

    delay = '%i' % (interval/10)
    call(['convert', '-loop', '0',
          '-layers', 'Optimize',
          '-delay', delay,
          'tmp/frame_*.png', args.output])

    # cleanup
    for i in range(len(phases)):
        os.remove('tmp/frame_%05i.png' % i)
    

@mlab.show
@mlab.animate(delay=interval)
def anim():
    phase = 0.0
    while True:
        phase += dphase

        dr = s*sin(phase)*args.amplitude
        x = (1.+dr)*sin(Th)*cos(Ph)
        y = (1.+dr)*sin(Th)*sin(Ph)
        z = (1.+dr)*cos(Th)

        m.mlab_source.set(x=x, y=y, z=z)
        yield

anim()
