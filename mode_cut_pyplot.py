#!/usr/bin/env python

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as pl
import numpy as np
from scipy.special import sph_harm, factorial
from tomso import adipls, io

fig = pl.figure(figsize=(6,6))
# ax = fig.add_subplot(111, projection='3d')  # equivalent?
ax = pl.gca(projection='3d')  # equivalent?
ell = 20
emm = 16
amax = 0.55    # np.sqrt((2.*ell+1.)/(4.*np.pi)*factorial(l-m)/factorial(l+m)) ?

Ntheta = 101
Nphi = 51

def get_colour(theta, phi):
    a = sph_harm(emm, ell, 
                 np.outer(phi, np.ones(len(theta))), 
                 np.outer(np.ones(len(phi)), theta)).real
    a = 0.5+0.5*a/amax
    return pl.cm.seismic(a)
    

# common keywords
kw = {'antialiased': True,
      # 'linewidth': 0.0,
      'cstride': 2,
      'rstride': 2}

# phi = np.linspace(0.5*np.pi, 2.0*np.pi, 100)
theta = np.linspace(0.0, np.pi, Ntheta)

# plotting as three separate segments avoids mplot3d getting confused about order
# x,y,z,a all have shape Nphi, Ntheta
phi = np.linspace(0.5*np.pi, 1.0*np.pi, Nphi)
x = np.outer(np.cos(phi), np.sin(theta))
y = np.outer(np.sin(phi), np.sin(theta))
z = np.outer(np.ones(np.size(phi)), np.cos(theta))
a = get_colour(theta, phi)
ax.plot_surface(x, y, z, facecolors=a, **kw)

phi = np.linspace(1.0*np.pi, 1.5*np.pi, Nphi)
x = np.outer(np.cos(phi), np.sin(theta))
y = np.outer(np.sin(phi), np.sin(theta))
z = np.outer(np.ones(np.size(phi)), np.cos(theta))
a = get_colour(theta, phi)
ax.plot_surface(x, y, z, facecolors=a, **kw)

phi = np.linspace(0.0*np.pi, 0.5*np.pi, Nphi)
x = np.outer(np.cos(phi), np.sin(theta))
y = np.outer(np.sin(phi), np.sin(theta))
z = np.outer(np.ones(np.size(phi)), np.cos(theta))
a = get_colour(theta, phi)
ax.plot_surface(x, y, z, facecolors=a, **kw)

fgong = io.load_fgong('data/modelS.fgong')
css, eigs = adipls.load_amde('data/modelS.amde')
I = np.where(css['ell']==ell)[0]
i = I[np.argmin((css['nu_Ri'][I]-3.)**2)]
r = eigs[i][:,0][::3]
y1 = eigs[i][:,1][::3]
rho = np.interp(r, fgong['var'][::-1,0]/fgong['glob'][1], fgong['var'][::-1,4])

# r = np.linspace(0.,1.,51)
# a = np.outer(np.sin(2.*np.pi*r), np.ones(len(theta)))
a = np.outer(r*rho**0.5*y1, np.ones(len(theta)))
a = a*sph_harm(emm, ell, 
               np.outer(np.ones(len(r)), np.ones(len(theta))), 
               np.outer(np.ones(len(r)), theta)).real
a = pl.cm.seismic(0.5+0.5*a/np.max(np.abs(a)))

x = np.outer(r, np.sin(theta))
y = 0.*x
z = np.outer(r, np.cos(theta))
ax.plot_surface(x, y, z, facecolors=a, **kw)

y = -np.outer(r, np.sin(theta))
x = 0.*y
z = np.outer(r, np.cos(theta))
ax.plot_surface(x, y, z, facecolors=a, **kw)

ax.set_xlim(-0.6,0.6)
ax.set_ylim(-0.6,0.6)
ax.set_zlim(-0.6,0.6)
pl.axis('off')
ax.view_init(0.,-80.)
pl.show()
