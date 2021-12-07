#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--figsize', type=float, nargs=2,
                    help="figure size, passed to rcParams['figure.figsize']")
parser.add_argument('--levels', type=int, default=20,
                    help="number of levels passed to contourf (default 100)")
parser.add_argument('--padding', type=float, default=0.01,
                    help="fractional padding between edge and circle (default=0.01)")
parser.add_argument('--cmap', type=str, default='viridis',
                    help="name of matplotlib colourmap to use (default='viridis')")
args = parser.parse_args()

if args.figsize:
    pl.rcParams['figure.figsize'] = args.figsize


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
    response.close()
    np.save('data/hmi_rot2d.npy', rot2d)

    response = urlopen('http://jsoc.stanford.edu/SUM86/D917240671/S00000/err.2d')
    err2d = np.loadtxt(response.readlines())
    response.close()
    np.save('data/hmi_err2d.npy', err2d)

    response = urlopen('http://jsoc.stanford.edu/SUM86/D917240671/S00000/rmesh.orig')
    rmesh = np.loadtxt(response.readlines())[::4]
    response.close()
    np.save('data/hmi_rmesh.npy', rmesh)

# rot2d has 49 columns, latitudes are 90-i*15/8; i starts at 0
lat = np.array([15./8.*i for i in np.arange(49)])/180.*np.pi
r, th = np.meshgrid(rmesh, lat)

ax = pl.subplot(111, projection='polar')
b = args.padding
pl.subplots_adjust(top=1-b, bottom=b, left=b, right=1-b)

data = rot2d.T[::-1]
data[err2d.T[::-1]/data>0.01] = np.nan
ax.contourf(th, r, data, args.levels, cmap=args.cmap)
ax.contourf(np.pi-th, r, data, args.levels, cmap=args.cmap)
ax.contourf(-th, r, data, args.levels, cmap=args.cmap)
ax.contourf(th-np.pi, r, data, args.levels, cmap=args.cmap)

# plot base of convection zone
th = np.linspace(0., 2.*np.pi, 401)
r = np.ones(len(th))*0.713
ax.plot(th, r, 'k--')

ax.set_xticklabels([])
ax.set_yticklabels([])
ax.grid(False)
pl.show()
