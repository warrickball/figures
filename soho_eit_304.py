#!/usr/bin/env python

import numpy as np
from scipy import ndimage
from matplotlib import pyplot as pl
from astropy.io import fits

# files found using VSO
url = 'https://sohodata.nascom.nasa.gov//archive/soho/private/data/processed/eit/lz/2016/11/efz20161101.131938'
datafile = 'data/soho_eit_304.fits'

try:
    data = fits.open(datafile)[0].data
except IOError:
    try:
        from urllib2 import urlopen
    except ImportError:
        from urllib.request import urlopen
        
    response = urlopen(url)
    with open(datafile, 'wb') as f:
        f.write(response.read())

    response.close()

    data = fits.open(datafile)[0].data

# tips for stripping everything but image:
# https://stackoverflow.com/a/9295367/1299112
data = data[2:-2, 2:-2]
size = (6.0, 6.0)
dpi = len(data)/size[0]
fig = pl.figure()
fig.set_size_inches(size)
ax = pl.Axes(fig, [0., 0., 1., 1.])
ax.set_axis_off()
fig.add_axes(ax)
pl.imshow(np.log10(data), vmin=2.917, vmax=3.0, cmap='hot')
pl.savefig('soho_eit_304.png', dpi=dpi)
pl.savefig('soho_eit_304.jpg', dpi=dpi)
pl.show()
