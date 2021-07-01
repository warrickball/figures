#!/usr/bin/env python3

from argparse import ArgumentParser
import numpy as np
import matplotlib.pyplot as pl
from scipy.optimize import curve_fit

parser = ArgumentParser()
parser.add_argument('--filters', type=str, default='UBVRI',
                    help="selection of photometric bands to plot "
                    "(default='UBVRI')")
parser.add_argument('--marker', type=str, default='.',
                    help="marker style (default='.')")
args = parser.parse_args()

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

if 'U' in args.filters:
    pl.errorbar(UB['phase'], multiharmonic(UB['phase'], *(A_V+A_BV))+UB['UB'], fmt=args.marker, yerr=UB['err'], label='U', c='C4');
if 'B' in args.filters:
    pl.errorbar(BV['phase'], multiharmonic(BV['phase'], *A_V)+BV['BV'], fmt=args.marker, yerr=BV['err'], label='B', c='C0');
if 'V' in args.filters:
    pl.errorbar(V['phase'], V['V'], fmt=args.marker, yerr=V['err'], label='V', c='C2');
if 'R' in args.filters:
    pl.errorbar(VR['phase'], multiharmonic(VR['phase'], *A_V)-VR['VR'], fmt=args.marker, yerr=VR['err'], label='R', c='C1');
if 'I' in args.filters:
    pl.errorbar(RI['phase'], multiharmonic(RI['phase'], *(A_V-A_VR))-RI['RI'], fmt=args.marker, yerr=RI['err'], label='I', c='C3');

pl.gca().invert_yaxis()
pl.xlabel('phase')
pl.suptitle('δ Cephei')
if len(args.filters) > 1:
    pl.legend()
    pl.ylabel('magnitude')
else:
    pl.ylabel('%s magnitude' % args.filters)

pl.show()
