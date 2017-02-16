#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from astropy.io import fits

# GOLF data
# Garcia R., Turck-Chieze S., Boumier P. et al. 2005, A&A, 442, 385 
# http://irfu.cea.fr/Phocea/file.php?class=astimg&file=1130/GOLF_velocity_series_mean_pm1_pm2.fits.gz

# VIRGO/SPM data
# Fr{\"o}hlich, C., Romero, J., Roth, H., et al.\ 1995, \solphys, 162, 101
# Frohlich, C., Andersen, B.~N., Appourchaux, T., et al.\ 1997, \solphys, 170, 1
# Jim{\'e}nez, A., Roca Cort{\'e}s, T., \& Jim{\'e}nez-Reyes, S.~J.\ 2002, \solphys, 209, 247
# http://irfu.cea.fr/Sap/Phocea/Vie_des_labos/Ast/ast_visu.php?id_ast=3581

longnames = ['GOLF_velocity_series_mean_pm1_pm2.fits.gz',
             'SPM_red_intensity_series.fits',
             'SPM_green_intensity_series.fits',
             'SPM_blue_intensity_series.fits']

shortnames = ['GOLF_mean_pm1_pm2.fits.gz',
              'SPM_red.fits',
              'SPM_green.fits',
              'SPM_blue.fits']

filenumbers = ['1130','3581','3581','3581']

data = {'GOLF': None, 'red': None, 'green': None, 'blue': None}
        
for key, shortname, longname, filenumber in zip(data.keys(), shortnames, longnames, filenumbers):
    try:
        data[key] = fits.open('data/%s' % shortname)[0].data
    except IOError:
        print('Downloading %s...' % shortname)
        import urllib2
        response = urllib2.urlopen('http://irfu.cea.fr/Phocea/file.php?class=astimg&file=%s/%s' % (filenumber, longname))
        with open('data/%s' % shortname, 'w') as f:
            f.write(response.read())
            
        data[key] = fits.open('data/%s' % shortname)[0].data

        
for key, value in data.iteritems():
    x = value
    dt = 60.0
    t = dt*np.arange(len(x))

    y = np.fft.rfft(x)
    df = 1./t[-1]
    f = np.arange(0.0, 0.5/dt+1.1*df, df)

    pl.loglog(1e3*f, np.abs(y)**2, lw=1.0)

pl.xlabel('frequency (mHz)')
pl.ylabel('power')
pl.show()
