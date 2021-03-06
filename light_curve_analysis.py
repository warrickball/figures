#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import fsolve
from astropy.timeseries import LombScargle

# get a smooth sawtooth with parametric curve
# (t,y) = (x-kcos(x)/2π, cos(2πx))
# where k=1 is maximum skewness
# so, given t, need to be able to solve for x

pl.rcParams['font.size'] = 14.0

TAU = 2.*np.pi

def func(x, t):
    return x-0.75*np.cos(TAU*x)/TAU-t

np.random.seed(123790186)

t = np.sort(np.random.uniform(low=0, high=3, size=20))
x = np.squeeze([fsolve(func, ti, args=(ti,)) for ti in t])
y = np.cos(TAU*x)

# y = np.sin(2.*np.pi*t)
f, p = LombScargle(t, y).autopower(samples_per_peak=20, nyquist_factor=2)
P = f[np.argmax(p)]

fig, axs = pl.subplots(nrows=3, ncols=1)

axs[1].plot(f, p, 'k');
for i in range(3):
    I = (i < t) & (t < i+1)
    axs[0].plot(t[I], y[I], 'o', label='cycle %i' % (i+1))
    axs[2].plot((t[I]/P)%1, y[I], 'o')

# axs[0].legend(['cycle 1', 'cycle 2', 'cycle 3'], loc='lower left');

a = [-0.06, 3.06, -1.1, 1.1]
axs[0].axis(a)
axs[0].set_xlabel('time')
axs[0].set_ylabel('magnitude')
axs[0].text(a[0]+0.5*(a[1]-a[0]), a[2] + 0.98*(a[3]-a[2]), 'light curve', va='top', ha='center')

a = [f.min()-0.06, f.max()+0.06, -0.05, 1.05]
axs[1].axis(a)
axs[1].set_xlabel('frequency (cycles/day)')
axs[1].set_ylabel('amplitude (mag)')
axs[1].text(a[0]+0.5*(a[1]-a[0]), a[2] + 0.98*(a[3]-a[2]), 'Lomb–Scargle periodogram', va='top', ha='center')

a = [-0.02, 1.02, -1.1, 1.1]
axs[2].axis(a)
axs[2].set_xlabel('phase')
axs[2].set_ylabel('magnitude')
axs[2].text(a[0]+0.5*(a[1]-a[0]), a[2] + 0.98*(a[3]-a[2]), 'phase-folded light curve', va='top', ha='center')

fig.tight_layout()
pl.show()
