#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
from numpy import sin, cos, pi
from argparse import ArgumentParser

parser = ArgumentParser(description="""Uses matplotlib to animates a spherical harmonic with a chosen
angular degree and azimuthal order.  """)
parser.add_argument('-l', '--ell', type=int, default=6,
                    help="angular degree")
parser.add_argument('-m', '--emm', type=int, default=3,
                    help="azimuthal order")
parser.add_argument('-o', '--output', type=str, default=None,
                    help="save figure to given filename without displaying "
                    "it (forces software rendering)")
parser.add_argument('--Ntheta', type=int, default=101,
                    help="number of points in latitude (default=101)")
parser.add_argument('--Nphi', type=int, default=101,
                    help="number of points in longitude (default=101)")
parser.add_argument('--cmap', type=str, default='seismic',
                    help="colour map for surface of sphere (default='seismic')")
parser.add_argument('--vmax', type=float, default=1.0,
                    help="maximum range of colour map; < 1.0 will saturate "
                    "(default=1.0)")
parser.add_argument('--vmin', type=float, default=None,
                    help="minimum range of colour map (default=-vmax)")
parser.add_argument('--pattern', type=str, default='',
                    help="pattern for surface colours; options are \n"
                    "'displacement' or 'dr' to use the local displacement "
                    "or anything else for a static spherical harmonic map.")
parser.add_argument('-a', '--amplitude', type=float, default=1.0/3.0,
                    help="amplitude of oscillation (default=1/3)")
parser.add_argument('-P', '--period', type=float, default=1.0,
                    help="period of oscillation, in seconds (default=1.0)")
parser.add_argument('--Nframes', type=int, default=40,
                    help="number of frames per oscillation (default=40)")
parser.add_argument('--view', type=float, nargs=2, default=[35.0, 45.0],
                    help="viewing angle")
args = parser.parse_args()

Nframes = args.Nframes
interval = args.period/Nframes*1e3  # in milliseconds
vmin = args.vmin if args.vmin else -args.vmax

def update(i, ax):
    ax.cla()

    phase = 2.*pi*i/Nframes
    dr = s*sin(phase)*args.amplitude
    x = (1.+dr)*sin(Th)*cos(Ph)
    y = (1.+dr)*sin(Th)*sin(Ph)
    z = (1.+dr)*cos(Th)
    v = s/(args.vmax-vmin)
    if args.pattern in ['displacement', 'dr']:
        v = v*np.sin(phase)

    v += 0.5

    surf = ax.plot_surface(x, y, z,
                           facecolors=pl.cm.get_cmap(args.cmap)(v),
                           **plot_kwargs)

    ax.set_xlim(-0.9,0.9)
    ax.set_ylim(-0.9,0.9)
    ax.set_zlim(-0.9,0.9)
    pl.axis('off')
    ax.view_init(*args.view)

    return surf,

plot_kwargs = {'rstride': 1,
               'cstride': 1,
               'linewidth': 0,
               'antialiased': False}

ell = args.ell
emm = args.emm

fig = pl.figure(figsize=(6,6))
# ax = Axes3D.Axes3D(fig)  # this is what tutorial uses
ax = pl.gca(projection='3d')
th = np.linspace(0., pi, args.Ntheta)
ph = np.linspace(-pi, pi, args.Nphi)
Th, Ph = np.meshgrid(th, ph)
s = sph_harm(emm,ell,Ph,Th).real
s = s/np.max(s)

update(0, ax)

ani = animation.FuncAnimation(fig, update, Nframes,
                              fargs=(ax,), interval=interval, repeat=True)

# Much smoother if we save it
if args.output:
    ani.save(args.output, writer='imagemagick')
else:
    pl.show()
