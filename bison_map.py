#!/usr/bin/env python

# based on example at
# http://scitools.org.uk/cartopy/docs/latest/examples/global_map.html

import matplotlib.pyplot as pl
import cartopy.crs as ccrs

ax = pl.axes(projection=ccrs.Robinson())
ax.set_global()

ax.stock_img()
ax.coastlines()

# co-ordinates taken from Wikipedia article on BiSON 2017-11-20 11:27
coords = {'Mount Wilson': (34.22, -118.07),
          'Las Campanas': (-29.01597, -70.69208),
          'Teide Observatory': (28.3, -16.5097),
          'SAAO': (-33.9347, 18.4776),
          'Carnarvon': (-24.869167, 113.704722),
          'Paul Wild Observatory': (-30.314, 149.562)}

for lat, lon in coords.values():
    pl.plot(lon, lat, 'o', transform=ccrs.PlateCarree())

pl.show()
