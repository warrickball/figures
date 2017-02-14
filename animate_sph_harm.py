#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm
from numpy import sin, cos, pi

period = 1.0  # in seconds
Nframes = 50
interval = period/Nframes*1e3  # is milliseconds

def update(i, ax):
    ax.cla()
    dr = sph_harm(emm,ell,Ph,Th).real*cos(2.*pi*i/Nframes)
    x = (1.+dr)*sin(Th)*cos(Ph)
    y = (1.+dr)*sin(Th)*sin(Ph)
    z = (1.+dr)*cos(Th)
    s = ax.plot_surface(x, y, z,
                        facecolors=pl.cm.seismic(0.5+dr),
                        **plot_kwargs)

    ax.set_xlim(-0.9,0.9)
    ax.set_ylim(-0.9,0.9)
    ax.set_zlim(-0.9,0.9)
    pl.axis('off')

    return s,

plot_kwargs = {'rstride':2,
               'cstride':2,
               'linewidth':0,
               'antialiased':True}

ell = 3
emm = 2

fig = pl.figure(figsize=(6,6))
# ax = Axes3D.Axes3D(fig)  # this is what tutorial uses
ax = pl.gca(projection='3d')

# dtheta = pi/100.
# dphi = 2.*pi/100.

# th = np.arange(0., pi+dtheta, dtheta)
# ph = np.arange(-pi, pi+dphi, dphi)
# Th, Ph = np.meshgrid(th, ph)
Th, Ph = np.meshgrid(np.linspace(0., pi, 51),
                     np.linspace(-pi, pi, 51))

update(0, ax)

ani = animation.FuncAnimation(fig, update, Nframes,
                              fargs=(ax,), interval=interval, repeat=False)

# Much smoother if we save it
# pl.show()
ani.save('animate_sph.gif', writer='imagemagick')
