#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl

# B. Mosser et al., bibcode 2014A&A...572L...5M

# similar data available in 
# B. Mosser et al., bibcode 2012A%26A...540A.143M
# ftp://cdsarc.u-strasbg.fr/pub/cats/J/A%2BA/540/A143

try:
    data = np.load('data/DPi_Dnu.npy')
except IOError:
    from astropy.io import ascii

    data = ascii.read('ftp://cdsarc.u-strasbg.fr/pub/cats/J/A%2BA/572/L5/table1.dat',
                      readme='ftp://cdsarc.u-strasbg.fr/pub/cats/J/A%2BA/572/L5/ReadMe').as_array()
    np.save('data/DPi_Dnu.npy', data)

data = data[np.argsort(data['Mass'])]
scat = pl.scatter(data['Dnu'], data['DPi1'], c=data['Mass'], s=15, lw=0,
                  cmap=pl.cm.jet)
ax = pl.gca()
ax.set_xscale('log')
ax.set_yscale('log')
pl.axis([2., 100., 40., 1000.])
# pl.xlabel(r'$\Delta\nu$')
# pl.ylabel(r'$\Delta\Pi_1$')
pl.xlabel('large separation ($\mu$Hz)')
pl.ylabel('period spacing (s)')
fig = pl.gcf()
cbax = fig.add_axes([0.14, 0.83, 0.5, 0.03])
pl.colorbar(scat, cax=cbax, orientation='horizontal')
pl.show()
