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
parser.add_argument('--figsize', type=float, nargs=2,
                    help="figure size, passed to rcParams['figure.figsize']")
parser.add_argument('--levels', type=int, default=100,
                    help="number of levels passed to contourf (default 100)")
parser.add_argument('--squeeze', action='store_const', default=False,
                    const=True, help="removes space around figure so that outer edge touches the figure border")
args = parser.parse_args()

if args.figsize:
    pl.rcParams['figure.figsize'] = args.figsize

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
    
r = eigs[i][:,0]
y1 = eigs[i][:,1]
rho = np.interp(r, var[::-1,0]/glob[1], var[::-1,4])
theta = np.linspace(0.0, 2.*np.pi, Ntheta)
a = np.outer(r*rho**0.5*y1, np.ones(len(theta)))
# use pi/2-theta instead of theta so that 0 degrees is at the top
a = a*sph_harm(args.emm, args.ell, 
               np.outer(np.ones(len(r)), np.ones(len(theta))), 
               np.outer(np.ones(len(r)), np.pi/2-theta)).real
# a = 0.5 + 0.5*a/np.max(np.abs(a))
amax = np.max(np.abs(a))

grid = np.meshgrid(theta, r)

ax = pl.subplot(111, projection='polar')
if args.squeeze:
    pl.subplots_adjust(top=1, bottom=0, left=0, right=1)

ax.contourf(grid[0], grid[1], a, args.levels, cmap='seismic',
            vmin=-amax, vmax=amax)
ax.set_xticklabels([])
ax.set_yticklabels([])
pl.show()
