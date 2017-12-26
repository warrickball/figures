#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from scipy.ndimage import imread
from argparse import ArgumentParser

parser = ArgumentParser(description="""
Plots a single HMI dopplegram, a difference between two consecutive
dopplergrams, or both (side-by-side).""")
parser.add_argument('option', help=" Choice of what to plot.  Options are `one`, `diff` or `both`.", type=str)
args = parser.parse_args()

def load(source, filename):
    try:
        data = imread(filename)
    except IOError:
        try:
            from urllib import urlretrieve
        except ImportError:
            from urllib.request import urlretrieve

        urlretrieve(source, filename)
        data = imread(filename)

    return data.astype('float')

def plot_one():
    data = np.copy(data1)
    N = len(data)
    x = np.linspace(-1., 1., N).reshape((1,-1))
    y = np.linspace(-1., 1., N).reshape((-1,1))
    z = x**2 + y**2
    data[(data < 15) & (z > 0.885)] = np.nan  # clip edges
    data[-100:, :N//2] = np.nan  # remove watermark
    pl.imshow(data, cmap='seismic')

def plot_diff():
    data = np.copy(data2-data1)
    N = len(data)
    x = np.linspace(-1., 1., N).reshape((1,-1))
    y = np.linspace(-1., 1., N).reshape((-1,1))
    z = x**2 + y**2
    data[(data < 4) & (z > 0.897)] = np.nan  # clip edges
    data[-100:, :N//2] = np.nan  # remove watermark
    pl.imshow(data, cmap='seismic')

source1 = 'https://sdo.gsfc.nasa.gov/assets/img/browse/2017/12/15/20171215_000438_4096_HMID.jpg'
filename1 = 'data/hmi_dopplergram1.jpg'
source2 = 'https://sdo.gsfc.nasa.gov/assets/img/browse/2017/12/15/20171215_002838_4096_HMID.jpg'
filename2 = 'data/hmi_dopplergram2.jpg'

if args.option == 'one':
    data1 = load(source1, filename1)
elif args.option == 'diff' or args.option == 'both':
    data1 = load(source1, filename1)
    data2 = load(source2, filename2)
else:
    raise ValueError('`option` must be `one`, `diff` or `both`')

# tips for stripping everything but image:
# https://stackoverflow.com/a/9295367/1299112
size = (6.0, 6.0)
try:
    dpi = len(data1)/size[0]
except NameError:
    dpi = len(data2)/size[0]
    
fig = pl.figure()

if args.option == 'both':
    fig.set_size_inches((2*size[0], size[1]))
    ax = pl.axes([0., 0., 0.5, 1.])
    ax.set_axis_off()
    plot_one()
    ax = pl.axes([0.5, 0., 0.5, 1.])
    ax.set_axis_off()
    plot_diff()
else:
    fig.set_size_inches(size)
    ax = pl.axes([0., 0., 1., 1.])
    ax.set_axis_off()
    if args.option == 'one':
        plot_one()
    else:
        plot_diff()

# pl.savefig('difference_dopplergram.jpg', dpi=dpi, quality=80)
pl.show()
