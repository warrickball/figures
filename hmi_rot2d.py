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

# adapted from
# http://matplotlib.org/examples/axes_grid/demo_floating_axes.html
fig = pl.gcf()
tr = PolarAxes.PolarTransform()
angle_ticks = [(0.0, r"$0^\circ$"),
               (0.25*np.pi, r"$45^\circ$"),
               (0.5*np.pi, r"$90^\circ$")]
grid_locator1 = FixedLocator([v for v, s in angle_ticks])
tick_formatter1 = DictFormatter(dict(angle_ticks))

grid_locator2 = MaxNLocator(2)

grid_helper = floating_axes.GridHelperCurveLinear(
    tr, extremes=(0, 0.5*np.pi, 0, 1),
    grid_locator1=grid_locator1,
    grid_locator2=grid_locator2,
    tick_formatter1=tick_formatter1,
    tick_formatter2=None)

ax1 = floating_axes.FloatingSubplot(fig, 111, grid_helper=grid_helper)
fig.add_subplot(ax1)

ax1.axis["left"].set_axis_direction("bottom")
ax1.axis["right"].set_axis_direction("top")

ax1.axis["bottom"].set_visible(False)
ax1.axis["top"].set_axis_direction("bottom")
ax1.axis["top"].toggle(ticklabels=True, label=True)
ax1.axis["top"].major_ticklabels.set_axis_direction("top")
ax1.axis["top"].label.set_axis_direction("top")

ax1.axis["left"].label.set_text(r'$r/R_\odot$')

aux_ax = ax1.get_aux_axes(tr)
aux_ax.patch = ax1.patch
ax1.patch.zorder = 0.9

Ncontours = 20
data = rot2d.T[::-1]
data[err2d.T[::-1]/data>0.01] = np.nan
c = aux_ax.contourf(th, r, data, Ncontours)
pl.colorbar(c, label='rotation rate (nHz)')
 
# plot base of convection zone
th = np.linspace(0., np.pi/2, 101)
r = np.ones(len(th))*0.713
aux_ax.plot(th, r, 'k--')

pl.show()
