#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl

# data from SDO/HMI webpage
# http://jsoc.stanford.edu/HMI/Global_products.html

try:
    data = np.load('data/hmi_modes.npy')
except IOError:
    import urllib2
    response = urllib2.urlopen('http://jsoc.stanford.edu/SUM99/D951224643/S00000/m10qr.8488')
    data = np.loadtxt(response.readlines())
    np.save('data/hmi_modes.npy', data)

l, n, obs, err = data[:,[0,1,2,7]].T
pl.semilogx(obs/(0.5+l), (n+1.5)/2./obs/1e-6, 'o')
# pl.scatter(obs/(0.5+l), (n+1.46)/2./obs/1e-6, s=10, c=n)
pl.xlabel(r'$\nu_{n\ell}/(\ell+1/2)\,(\mu\mathrm{Hz})$')
pl.ylabel(r'$(n+\alpha)/2\nu_{n\ell}\,(\mathrm{s})$')
pl.show()
