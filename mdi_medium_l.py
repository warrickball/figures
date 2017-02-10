#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from astropy.io import fits

# E.J. Rhodes, Jr., A.G. Kosovichev, P.H. Scherrer, J. Schou & J. Reiter
# https://soho.nascom.nasa.gov/publications/CDROM1/papers/rhodes/lnu0.fts

try:
    data = fits.open('data/mdi_medium_l.fits')[0].data
except IOError:
    import urllib2
    response = urllib2.urlopen('https://soho.nascom.nasa.gov/publications/CDROM1/papers/rhodes/lnu0.fts')
    with open('data/mdi_medium_l.fits','w') as f:
        f.write(response.read())
        
    data = fits.open('data/mdi_medium_l.fits')[0].data
    
f_max = 0.0803755*103680/1e6  # according to website
kwargs = {'aspect':'auto',
          'origin':'lower',
          'interpolation':'none',
          'vmax':np.max(np.log10(data))-np.log10(100.0),
          'extent':[0.,300.,0.,f_max*1e3]}

pl.imshow(np.log10(data), **kwargs)
pl.xlabel(r'angular degree $\ell$')
pl.ylabel(r'frequency (mHz)')
pl.show()

