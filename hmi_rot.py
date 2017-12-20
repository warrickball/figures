#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl

from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist import floating_axes
from mpl_toolkits.axisartist.grid_finder import (FixedLocator, MaxNLocator,
                                                 DictFormatter)
# data from SDO/HMI webpage
# http://jsoc.stanford.edu/HMI/Global_products.html

try:
    rot2d = np.load('data/hmi_rot2d.npy')
    err2d = np.load('data/hmi_err2d.npy')
    rmesh = np.load('data/hmi_rmesh.npy')
except IOError:
    try:
        from urllib2 import urlopen
    except ImportError:
        from urllib.request import urlopen
        
    response = urlopen('http://jsoc.stanford.edu/SUM86/D917240671/S00000/rot.2d')
    rot2d = np.loadtxt(response.readlines())
    np.save('data/hmi_rot2d.npy', rot2d)

    response = urlopen('http://jsoc.stanford.edu/SUM86/D917240671/S00000/err.2d')
    err2d = np.loadtxt(response.readlines())
    np.save('data/hmi_err2d.npy', err2d)

    response = urlopen('http://jsoc.stanford.edu/SUM86/D917240671/S00000/rmesh.orig')
    rmesh = np.loadtxt(response.readlines())[::4]
    np.save('data/hmi_rmesh.npy', rmesh)

incs = np.array([90.0-i*15./8. for i in range(len(rot2d[0]))])
# rot2d has 49 columns, latitudes are 90-i*15/8; i starts at 0
for inc, rot, err in list(zip(incs, rot2d.T, err2d.T))[8::8][::-1]:
    c, = pl.plot(rmesh, rot, label='$%i^\circ$' % inc)
    pl.fill_between(rmesh, rot-err, rot+err,
                    color=c.get_color(), alpha=0.5)

# pl.plot([0.713,0.713], [np.min(rot2d), np.max(rot2d)], 'k--')
pl.xlabel(r'$r/R_\odot$')
pl.ylabel(r'$\Omega/\mathrm{nHz}$')
pl.legend()
pl.show()
