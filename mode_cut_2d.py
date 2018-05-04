#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
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
args = parser.parse_args('')

if args.freq and args.enn:
        raise ValueError("""you can only use one of the radial order (-n,
        --enn) or the cyclic frequency (-f, --freq)""")


amax = 0.55    # np.sqrt((2.*ell+1.)/(4.*np.pi)*factorial(l-m)/factorial(l+m)) ?

Ntheta = 201

def get_colour(theta, phi):
    a = sph_harm(args.emm, args.ell, 
                 np.outer(phi, np.ones(len(theta))), 
                 np.outer(np.ones(len(phi)), theta)).real
    a = 0.5+0.5*a/amax
    return pl.cm.seismic(a)
    
glob, var = fgong.load_fgong('data/modelS.fgong')
css, eigs = adipls.load_amde('data/modelS.amde')
I = np.where(css['ell']==args.ell)[0]
if args.freq:
    i = I[np.argmin((css['nu_Ri'][I]-args.freq)**2)]
elif args.enn:
    i = I[css['enn'][I]==args.enn][0]
else:
    i = I[css['enn'][I]==14][0]
    
r = eigs[i][:,0][::3]
y1 = eigs[i][:,1][::3]
rho = np.interp(r, var[::-1,0]/glob[1], var[::-1,4])
theta = np.linspace(0.0, 2.*np.pi, Ntheta)
a = np.outer(r*rho**0.5*y1, np.ones(len(theta)))
a = a*sph_harm(args.emm, args.ell, 
               np.outer(np.ones(len(r)), np.ones(len(theta))), 
               np.outer(np.ones(len(r)), theta)).real
# a = 0.5 + 0.5*a/np.max(np.abs(a))
amax = np.max(np.abs(a))

grid = np.meshgrid(theta, r)

ax = pl.subplot(111, projection='polar')
ax.contourf(grid[0], grid[1], a, 100, cmap='seismic',
            vmin=-amax, vmax=amax)
pl.show()
