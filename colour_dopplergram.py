#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from scipy.ndimage import imread

source = 'https://sdo.gsfc.nasa.gov/assets/img/browse/2017/12/15/20171215_000438_4096_HMID.jpg'
filename = 'data/hmi_dopplergram.jpg'

try:
    data = imread(filename)
except IOError:
    try:
        from urllib import urlretrieve
    except ImportError:
        from urllib.request import urlretrieve

    urlretrieve(source, filename)
    data = imread(filename)

# tips for stripping everything but image:
# https://stackoverflow.com/a/9295367/1299112
size = (6.0, 6.0)
dpi = len(data)/size[0]
fig = pl.figure()
fig.set_size_inches(size)
ax = pl.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)

N = len(data)
data = np.array(data, dtype='float')
x = np.linspace(-1., 1., N).reshape((1,-1))
y = np.linspace(-1., 1., N).reshape((-1,1))
z = x**2 + y**2
data[(data < 15) & (z > 0.88)] = np.nan  # clip edges
data[-100:, :N//2] = np.nan  # remove watermark
pl.imshow(data, cmap='seismic')
# pl.savefig('colour_dopplergram.jpg', dpi=dpi, quality=80)
pl.show()
