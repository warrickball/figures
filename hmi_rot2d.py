#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl

# data from SDO/HMI webpage
# http://jsoc.stanford.edu/HMI/Global_products.html

try:
    rot2d = np.load('data/hmi_rot2d.npy')
    err2d = np.load('data/hmi_err2d.npy')
    rmesh = np.load('data/hmi_rmesh.npy')
except IOError:
    import urllib2
    response = urllib2.urlopen('http://jsoc.stanford.edu/SUM86/D917240671/S00000/rot.2d')
    rot2d = np.loadtxt(response.readlines())
    np.save('data/hmi_rot2d.npy', rot2d)

    response = urllib2.urlopen('http://jsoc.stanford.edu/SUM86/D917240671/S00000/err.2d')
    err2d = np.loadtxt(response.readlines())
    np.save('data/hmi_err2d.npy', err2d)

    response = urllib2.urlopen('http://jsoc.stanford.edu/SUM86/D917240671/S00000/rmesh.orig')
    rmesh = np.loadtxt(response.readlines())[::4]
    np.save('data/hmi_rmesh.npy', rmesh)

# rot2d has 49 columns, latitudes are 90-i*15/8; i starts at 0
lat = np.array([15./8.*i for i in np.arange(49)])/180.*np.pi
r, th = np.meshgrid(rmesh, lat)

fig = pl.figure()
ax = pl.subplot(111, polar=True)

# I currently fill out the full circle
# probably better to just do first quadrant (so symmetry is clear)
# but this happens to be hard in matplotlib...
# http://stackoverflow.com/questions/27433310/how-to-clip-polar-plot-in-pylab-pyplot
Ncontours = 100
pl.contourf(th, r, rot2d.T, Ncontours)
pl.contourf(np.pi-th, r, rot2d.T, Ncontours)
pl.contourf(np.pi+th, r, rot2d.T, Ncontours)
pl.contourf(2.*np.pi-th, r, rot2d.T, Ncontours)

pl.show()
