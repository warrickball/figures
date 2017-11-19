#!/usr/bin/env python

"""Intensity and velocity power spectra of the Sun, showing the
low-degree oscillations.  Data are taken from GOLF and VIRGO/SPM.

"""


import numpy as np
from matplotlib import pyplot as pl
from astropy.io import fits
import argparse

# GOLF data
# Garcia R., Turck-Chieze S., Boumier P. et al. 2005, A&A, 442, 385 
# http://irfu.cea.fr/Phocea/file.php?class=astimg&file=1130/GOLF_velocity_series_mean_pm1_pm2.fits.gz

# VIRGO/SPM data
# Fr{\"o}hlich, C., Romero, J., Roth, H., et al.\ 1995, \solphys, 162, 101
# Frohlich, C., Andersen, B.~N., Appourchaux, T., et al.\ 1997, \solphys, 170, 1
# Jim{\'e}nez, A., Roca Cort{\'e}s, T., \& Jim{\'e}nez-Reyes, S.~J.\ 2002, \solphys, 209, 247
# http://irfu.cea.fr/Sap/Phocea/Vie_des_labos/Ast/ast_visu.php?id_ast=3581

parser = argparse.ArgumentParser()
parser.add_argument('datasets', type=str, nargs='+',
                    help="List of datasets to be plotted.  Options are " + \
                    "'golf', 'spm_red,' 'spm_green' or 'spm_blue'.")
parser.add_argument('-s','--smoothing', type=int, nargs=1, default=20,
                    help="Number of points for boxcar smoothing.")
parser.add_argument('--legend', action='store_const', const=True,
                    default=False, help="Add a legend to the plot.")
args = parser.parse_args()

# the Fourier transforms take some time to calculate, so I do some
# bookkeepping to make sure we only compute them once
fftnames = {'golf':      'GOLF_fft.npy',
            'spm_red':   'SPM_red_fft.npy',
            'spm_green': 'SPM_green_fft.npy',
            'spm_blue':  'SPM_blue_fft.npy'}
longnames = {'golf':      'GOLF_velocity_series_mean_pm1_pm2.fits.gz',
             'spm_red':   'SPM_red_intensity_series.fits',
             'spm_green': 'SPM_green_intensity_series.fits',
             'spm_blue':  'SPM_blue_intensity_series.fits'}
shortnames = {'golf':      'GOLF_mean_pm1_pm2.fits.gz',
              'spm_red':   'SPM_red.fits',
              'spm_green': 'SPM_green.fits',
              'spm_blue':  'SPM_blue.fits'}
filenumbers = {'golf':      '1130',
               'spm_red':   '3581',
               'spm_green': '3581',
               'spm_blue':  '3581'}
names = {'golf':      'GOLF',
         'spm_red':   'SPM red',
         'spm_green': 'SPM green',
         'spm_blue':  'SPM blue'}
colors = {'golf':      'k',
          'spm_red':   'r',
          'spm_green': 'g',
          'spm_blue':  'b'}

dt = 60.0  # cadence, which happens to be the same for all timeseries
for key in args.datasets:
    k = key.lower()
    try:
        y = np.load('data/%s' % fftnames[k])
    except IOError:
        try:
            timeseries = fits.open('data/%s' % shortnames[k])[0].data
        except IOError:
            print('Downloading %s...' % shortnames[k])
            import urllib2
            response = urllib2.urlopen('http://irfu.cea.fr/Phocea/file.php?class=astimg&file=%s/%s' % (filenumber[k], longname[k]))
            with open('data/%s' % shortnames[k], 'w') as f:
                f.write(response.read())
            
            timeseries = fits.open('data/%s' % shortnames[k])[0].data

        y = np.fft.rfft(timeseries)
        np.save('data/%s' % fftnames[k], y)

    t = dt*np.arange(2*len(y)-2)
    df = 1./t[-1]
    f = np.arange(0.0, 0.5/dt+1.1*df, df)
    yy = np.convolve(np.abs(y)**2, np.ones(args.smoothing)/args.smoothing, mode='same')/np.median(np.abs(y[np.abs(f)-7.8<0.2])**2)

    pl.loglog(1e3*f, yy, lw=1.0, label=names[k], color=colors[k])

    if args.legend: pl.legend()

pl.xlabel('frequency (mHz)')
pl.ylabel('power relative to background')
pl.show()
