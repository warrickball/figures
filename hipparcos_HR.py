#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--annotate', action='store_true')
args = parser.parse_args()

from matplotlib import pyplot as pl
import numpy as np
from scipy.interpolate import interp1d
            
try:
    data = np.load('data/hip_main.npy')
except IOError:
    from astroquery.vizier import Vizier

    data = Vizier(columns=['**', '+_r'], row_limit=-1).get_catalogs('I/311/hip2')[0].as_array().data
    np.save('data/hip_main.npy', data)

data = data[data['e_Hpmag']<0.1]
data = data[data['Plx']>0]
data = data[data['e_Plx']/data['Plx']<0.05]

Hp = data['Hpmag'] + 5.0 + 5.0*np.log10(1e-3*data['Plx'])
BV = data['B-V']

I = ~np.isnan(BV)
BV = BV[I]
Hp = Hp[I]
# pl.hexbin(BV, Hp, cmap='gray_r', gridsize=50, bins='log')
# pl.hist2d(BV, Hp, cmap='gray_r', bins=100)
pl.plot(BV, Hp, 'ko', ms=4, alpha=0.1, mew=0)

# original annotations
if args.annotate:
    MS_BV, MS_Hp = np.array([[-0.18,-2.46], [0.0,0.93], [0.2,2.45],
                             [0.4,3.3], [0.63,5.], 
                             [0.9,6.4], [1.2,7.8],
                             [1.37,8.48], [1.57,11.75]]).T
    MS_BVi = np.linspace(np.min(MS_BV), np.max(MS_BV), 100)
    pl.plot(MS_BVi, interp1d(MS_BV, MS_Hp, kind='cubic')(MS_BVi), 'g-', lw=10,
            alpha=0.6, solid_capstyle='round')  # MS
    pl.fill_between([0.25,0.75,0.85],[0.,2.8,2.0], [0.,0.,2.], alpha=0.25,
                    lw=0)  # Hertzsprung gap
    pl.plot([0.85,1.6], [3.6,-1.05], 'r-', lw=25, alpha=0.2,
            solid_capstyle='round')  # RGB
    pl.plot([0.92,1.07],[0.9,1.],'r-', lw=25, alpha=0.2,
            solid_capstyle='round')  # RC
    pl.plot([-0.5,1.0], [8.6,17.34], 'b-', lw=25, alpha=0.1)  # WD
    pl.plot(0.653, 4.86, 'y*', ms=15)  # Sun

pl.axis([-0.5, 2., 15., -5.])
pl.xlabel(r'$B-V$')
pl.ylabel(r'$M_\mathrm{Hipparcos}$')
pl.show()
