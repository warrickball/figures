#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit

P = 5.366249    # Wikipedia
E = 2455479.905 # Paper

def JD_to_phase(JD, P, E):
    return ((JD-E)/P + 0.5)%1 - 0.5

def multiharmonic(t, *A):
    n = (len(A)-1)//2
    y = A[0]*np.ones_like(t)
    for i in range(n):
        y += A[i+1]*np.sin(2.*np.pi*(i+1)*t) + A[n+i+1]*np.cos(2.*np.pi*(i+1)*t)
        
    return y

def build_fit(t, y, n=10):
    r = [[np.nanmean(y)]]
    for i in range(1, 10):
        r = curve_fit(multiharmonic, t, y, np.hstack([r[0], 0., 0.]))

    return r[0]

# data extracted manually from
# Engle et al. 2014, ApJ, 794, 80
# https://ui.adsabs.harvard.edu/abs/2014ApJ...794...80E
V, UB, BV, VR, RI = [np.genfromtxt('data/del_cep_%s.dat' % name, names=True, skip_header=3)
                     for name in ['V', 'U-B', 'B-V', 'V-R', 'R-I']]

# P and E don't differ enough to matter for plot
# r = curve_fit(JD_to_phase,
#               np.hstack([d['HJD'] for d in [UB, BV, V, VR, RI]]),
#               np.hstack([d['phase'] for d in [UB, BV, V, VR, RI]]),
#               p0=(P,E), full_output=True)
# P, E = r[0]

A_V = build_fit(V['phase'], V['V'])
A_BV = build_fit(BV['phase'], BV['BV'])
A_UB = build_fit(UB['phase'], UB['UB'])
A_VR = build_fit(VR['phase'], VR['VR'])
A_RI = build_fit(RI['phase'], RI['RI'])

pl.errorbar(UB['phase'], multiharmonic(UB['phase'], *(A_V+A_BV))+UB['UB'], fmt='.', yerr=UB['err'], label='U', c='C4');
pl.errorbar(BV['phase'], multiharmonic(BV['phase'], *A_V)+BV['BV'], fmt='.', yerr=BV['err'], label='B', c='C0');
pl.errorbar(V['phase'], V['V'], fmt='.', yerr=V['err'], label='V', c='C2');
pl.errorbar(VR['phase'], multiharmonic(VR['phase'], *A_V)-VR['VR'], fmt='.', yerr=VR['err'], label='R', c='C1');
pl.errorbar(RI['phase'], multiharmonic(RI['phase'], *(A_V-A_VR))-RI['RI'], fmt='.', yerr=RI['err'], label='I', c='C3');
pl.gca().invert_yaxis()
pl.xlabel('phase')
pl.ylabel('magnitude')
pl.suptitle('δ Cephei')
pl.legend()
pl.show()
