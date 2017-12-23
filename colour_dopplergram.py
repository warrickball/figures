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

# data[data==0] = 120
    
# tips for stripping everything but image:
# https://stackoverflow.com/a/9295367/1299112
size = (6.0, 6.0)
dpi = len(data)/size[0]
fig = pl.figure()
fig.set_size_inches(size)
ax = pl.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)

pl.imshow(data, cmap='seismic')
pl.show()
