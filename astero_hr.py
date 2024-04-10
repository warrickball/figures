#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import UnivariateSpline
from scipy.optimize import root_scalar

from matplotlib.ticker import ScalarFormatter
from matplotlib.path import Path
import matplotlib.patches as patches

import mistery

# pl.rcParams['axes.formatter.min_exponent'] = 6

def polygon(xy, steps=5, **kwargs):
    # shape must be convex
    ax = pl.gca()
    midx, midy = np.mean(xy, axis=0)
    pts = [(x*(1-t)+t*midx, y*(1-t)+t*midy) for t in np.arange(0, 1, 1/steps) for x, y in xy]
    return patches.Polygon(pts, joinstyle='round', fill=False, **kwargs)

def add(x, y, *args, **kwargs):
    ax.loglog(10**np.squeeze(x), 10**np.squeeze(y), *args, **kwargs)

try:
    WD = np.load('data/astero_hr_WD_track.npy')
except:
    print("Downloading track for 0.6 Msun WD...")
    h = mistery.get_track(M=2.50)
    h = h[h['phase']>=0]
    WD = h[h['phase']==6]
    np.save('data/astero_hr_WD_track.npy', WD)

ZAMS_data = []
TAMS_data = []
for mass100 in [30, 60, 100, 200, 300, 600, 1000, 2000, 3000]: # , 6000, 10000]:
    filename = 'data/astero_hr_track_M%05i.npy' % mass100
    try:
        h = np.load(filename)
    except:
        print("Downloading track for %.2f Msun..." % (mass100*0.01))
        h = mistery.get_track(M=mass100*0.01)
        np.save(filename, h)

    ZAMS_data.append(h[h['phase']==0][0])
    TAMS_data.append(h[h['phase']==0][-1])

ZAMS_data = np.hstack(ZAMS_data)
TAMS_data = np.hstack(TAMS_data)
ZAMS = UnivariateSpline(ZAMS_data['log_L'], ZAMS_data['log_Teff'])
TAMS = UnivariateSpline(TAMS_data['log_L'], TAMS_data['log_Teff'])

highlight = dict(lw=20, alpha=0.5, zorder=-10, solid_capstyle='round')

fig, ax = pl.subplots()

# data

logL = np.linspace(-1, 5, 121)
# logL = np.linspace(ZAMS_data['log_L'].min(), ZAMS_data['log_L'].max(), 101)
ax.loglog(10**ZAMS(logL), 10**logL, 'k-')
# logL = np.linspace(TAMS_data['log_L'].min(), TAMS_data['log_L'].max(), 101)
# ax.plot(10**TAMS(logL), 10**logL, 'k-')

# main sequence

# use dlogL/dlogTeff = -20 for instability strip, adjust based on blue loops
# find where this intersects classical instability region on ZAMS (basically δ Scts)
# e.g. 7000--9000 K ⇒ logT = 3.85--3.95
# if logL₀ corresponds to some temperature logT₀,
# line is logL=20(logT-logT₀) + logL₀
logLr = root_scalar(lambda x: ZAMS(x)-3.85, x0=0.4, x1=0.6).root
ax.loglog([10**3.85, 10**(3.85+(logLr-5)/20)], [10**logLr, 10**5], 'k--')
logLb = root_scalar(lambda x: ZAMS(x)-3.95, x0=0.8, x1=1.0).root
ax.loglog([10**3.95, 10**(3.95+(logLb-5)/20)], [10**logLb, 10**5], 'k--')

# WD cooling track
ax.loglog(10**WD['log_Teff'], 10**WD['log_L'], 'k-')
p = np.polyfit(WD['log_Teff'][-10:], WD['log_L'][-10:], 1)
logT = np.linspace(3.8, WD['log_Teff'][-10])
ax.loglog(10**logT, 10**np.polyval(p, logT), 'k-')

# tracks
for mass100 in [100, 200, 300, 600, 1000]: # , 6000, 10000]:
    filename = 'data/astero_hr_track_M%05i.npy' % mass100
    try:
        h = np.load(filename)
    except:
        print("Downloading track for %.2f Msun..." % (mass100*0.01))
        h = mistery.get_track(M=mass100*0.01)
        np.save(filename, h)

    h = h[(0<=h['phase'])&(h['phase']<5)]
    ax.loglog(10**h['log_Teff'], 10**h['log_L'], 'k-', lw=1)

# highlights

# β Cep and SPB
logL = np.linspace(3.2, 4.8, 151)
add(ZAMS(logL), logL, **highlight)
logL = np.linspace(2.0, 3.5, 151)
add(ZAMS(logL), logL, **highlight)

# δ Sct, Cepheids and RR Lyrae
logL = np.linspace(logLr+0.2, logLb-0.2, 11)
add(ZAMS(logL), logL, **highlight)
f = lambda x: (logLr+logLb)/2 - 20*(x-3.9)
add([3.8, 3.7], [f(3.8), f(3.7)], **highlight);
add([3.8499, 3.8501], [f(3.8499), f(3.8501)], **highlight);

# solar-like, semi-regulars, Miras
# go parallel to instability strip but redder
f = lambda x: -0.1 - 20*(x-ZAMS(-0.1)) # sets where RG, SR & Mira regions intersect ZAMS
# add([ZAMS(0.25), ZAMS(-0.4), ZAMS(-0.1), 3.67], [0.25, -0.4, -0.1, f(3.67)], **highlight)
add([ZAMS(-0.4), ZAMS(0.25), 3.74, 3.67], [-0.4, 0.25, f(3.74), f(3.67)], **highlight)
add([3.66, 3.6], [f(3.66), f(3.6)], **highlight)
add([3.59, 3.5], [f(3.59), f(3.5)], **highlight)
# ax.add_patch(polygon([(3.59, f(3.59)), (3.5, f(3.5))

# WD instability regions
DOV = WD[WD['log_Teff']>5.1]
add(DOV['log_Teff'], DOV['log_L'], **highlight)
add([4.4, 4.6], np.polyval(p, [4.4, 4.6]), **highlight)
add([3.95, 4.15], np.polyval(p, [3.95, 4.15]), **highlight)

# subdwarf variables
add([4.3, 4.4], [1.4, 1.4], **highlight)
add([4.43, 4.5], [1.4, 1.4], **highlight)

# crude solar symbol
ax.scatter(5772., 1.0, fc='w', ec='k', zorder=10);
ax.scatter(5772., 1.0, s=1, fc='w', ec='k', zorder=10);

# ax.invert_xaxis()
ax.xaxis.set_major_formatter(ScalarFormatter())
ax.set_xticks([5e3, 1e4, 2e4, 5e4, 1e5, 2e5])
ax.set_yticks(10.0**np.arange(-5,7,2), minor=True, labels=[])
ax.set_xlabel('effective temperature (kelvin)')
ax.set_ylabel('luminosity (solar units)')
ax.tick_params(which='both', top=True, right=True)
ax.set_xlim(195000, 2500)
ax.set_ylim(1e-5, 3e5)

pl.show()

raise ValueError
# unrecognisable cruft

try:
    iso = np.load('data/astero_hr_iso.npy')
except:
    print("Downloading isochrone...")
    iso = mistery.get_isochrone(t=8.0)
    # iso = mistery.get_isochrone(t=8.0, FeH=-4)
    np.save('data/astero_hr_iso.npy', iso)

ax.plot(iso['log_Teff'], iso['log_L'])

# ZAMS from isochrones somewhy limited to max mass ~13 Msun
ZAMS = []
for t in np.unique(isos['log10_isochrone_age_yr']):
    iso = isos[(isos['log10_isochrone_age_yr']==t)&(isos['phase']==-1)]

    if len(iso) > 0:
        ZAMS.append(iso[-1])

ZAMS = np.hstack(ZAMS)

names = ['solar-like-ms', 'solar-like-rg', 'gdor', 'dsct',
         'spb', 'bcep', 'sdBVp', 'sdBVg', 'DOV', 'DBV', 'DAV',
         'rrlyr', 'dcep', 'mira']

# raw numbers from Aerts, JCD & Kurtz 2010
logT = np.vstack([[3.70, 3.70, 3.83, 3.82, 4.05, 4.25, 4.2, 4.4, 4.8, 4.4, 3.95, 3.88, 3.85, 3.75],
                  [3.82, 3.65, 3.90, 3.95, 4.35, 4.50, 4.5, 4.6, 5.1, 4.6, 4.15, 3.78, 3.55, 3.45]])
logL = np.vstack([[-0.5, -0.5, 0.7, 0.6, 2.0, 3.2, 1.2, 1.2, 1.5, -1.0, -2.6, 1.4, 2.0, 2.5],
                  [ 1.0,  2.0, 1.1, 2.0, 4.0, 5.0, 2.2, 2.6, 3.5,  0.7, -2.2, 1.7, 5.5, 4.0]])

ax.plot(logT, logL, 'o-');
pl.legend(names);



ax.plot(10**ZAMS['log_Teff'], ZAMS['log_L'])
ax.plot(10**TAMS['log_Teff'], TAMS['log_L'])
ax.plot(10**WD['log_Teff'], WD['log_L'])
p = np.polyfit(WD['log_Teff'][-10:], WD['log_L'][-10:], 1)
T = np.linspace(4, WD['log_Teff'][-10])
ax.plot(10**T, np.polyval(p, T))

# trying to get an instability strip
# logR = 0.680(17)logP + 1.146(25)
# Mv ~ 4.74 - 2.5logL = -2.43(12)(logP-1)-4.05(2)
# logL = (4.74 + 4.05(2) + 2.43(12)((logR+1.146(25))/0.680(17)-1))/2.5
#      = 4.1821 + 0.457*logR
#      = 4.1821 + 0.2287logL - 0.9148logT
# logL = 5.4221 - 1.186logT

# Sandage 2008: logT = 3.922-0.054logL
# reading off graph in ACDK, edges are (not quite parallel):
# ax.plot([3.65, 3.85], [5., 0.7], 'k-'); # red, gradient = -21.5
# ax.plot([3.7, 3.92], [5.25, 1.0], 'k-'); # blue, gradient =  -19.31818...
# intersect MS at ~7000 K and 8500K (logT ~ 3.85, 3.93)?
