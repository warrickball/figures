#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from astropy.io import fits

# Garcia R., Turck-Chieze S., Boumier P. et al. 2005, A&A, 442, 385 
# http://irfu.cea.fr/Phocea/file.php?class=astimg&file=1130/GOLF_velocity_series_mean_pm1_pm2.fits.gz
try:
    data = fits.open('data/golf_mean_pm1_pm2.fits.gz')[0].data
except IOError:
    import urllib2
    response = urllib2.urlopen('http://irfu.cea.fr/Phocea/file.php?class=astimg&file=1130/GOLF_velocity_series_mean_pm1_pm2.fits.gz')
    with open('data/golf_mean_pm1_pm2.fits.gz', 'w') as f:
        f.write(response.read())

    data = fits.open('data/golf_mean_pm1_pm2.fits.gz')[0].data
        
x = data
dt = 60.0
t = dt*np.arange(len(x))

y = np.fft.rfft(x)
df = 1./t[-1]
f = np.arange(0.0, 0.5/dt+1.1*df, df)

pl.loglog(1e3*f, np.abs(y)**2, lw=1.0)
pl.xlabel('frequency (mHz)')
pl.ylabel('power')
pl.show()
