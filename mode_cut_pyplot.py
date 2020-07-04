#!/usr/bin/env python

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as pl
import numpy as np
from scipy.special import sph_harm, factorial
from tomso import adipls, fgong
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-l', '--ell', type=int, default=20,
                    help="Angular degree l (default l=20)")
parser.add_argument('-m', '--emm', type=int, default=16,
                    help="Azimuthal order m (default m=16)")
parser.add_argument('-n', '--enn', type=int, default=None,
                    help="""Radial order n.  You can only choose one of this or the cyclic
                    frequency n (-f, --freq).  (default n=14)""")
parser.add_argument('-f', '--freq', type=float, default=None,
                    help="""Cyclic frequency in mHz.  You can only choose one of this or the
                    radial order n (-n, --enn).""")
args = parser.parse_args()

if args.freq and args.enn:
        raise ValueError("""you can only use one of the radial order (-n,
        --enn) or the cyclic frequency (-f, --freq)""")

fig = pl.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')  # equivalent?
# ax = pl.gca(projection='3d')  # equivalent?
fig.subplots_adjust(left=0., right=1., bottom=0., top=1.)
amax = 0.55    # np.sqrt((2.*ell+1.)/(4.*np.pi)*factorial(l-m)/factorial(l+m)) ?

Ntheta = 201
Nphi = 101

def get_colour(theta, phi):
    a = sph_harm(args.emm, args.ell, 
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

S = fgong.load_fgong('data/modelS.fgong', return_object=True)
amde = adipls.load_amde('data/modelS.amde', return_object=True)
I = np.where(amde.l==args.ell)[0]
if args.freq:
    i = I[np.argmin((amde.nu_Ri[I]-args.freq/1e3)**2)]
elif args.enn:
    i = I[amde.n[I]==args.enn][0]
else:
    i = I[amde.n[I]==14][0]
    
print('   n = %i' % amde.n[i])
print('freq = %.6f mHz' % (amde.nu_Ri[i]*1e3))
r = amde.eigs[i][:,0][::3]
y1 = amde.eigs[i][:,1][::3]
rho = np.interp(r, S.x[::-1], S.rho[::-1])

# r = np.linspace(0.,1.,51)
# a = np.outer(np.sin(2.*np.pi*r), np.ones(len(theta)))
a = np.outer(r*rho**0.5*y1, np.ones(len(theta)))
a = a*sph_harm(args.emm, args.ell, 
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
ax.view_init(0.,-75.)
pl.show()
