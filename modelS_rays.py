#!/usr/bin/env python

"""Illustrates (approximate) ray propagation for p-mode in Model S.
Currently only seems to work up to l=49 (and fails on a few
intermediate values too).  Still need to
- check that the inner turning point is correct,
- start integration from upper reflecting boundary, and
- add options for how many times to "bounce" at the surface.

"""

import numpy as np
from matplotlib import pyplot as pl
from tomso import io
from scipy import integrate as spint
from argparse import ArgumentParser

pl.rcParams['figure.figsize'] = 6, 6

parser = ArgumentParser()
parser.add_argument('--ell', type=int, nargs='+', default=[2, 20, 25, 75],
                    help='angular degree(s) of the desired rays')
parser.add_argument('--nu', type=float, default=3.002e-3,
                    help='cyclic frequency in Hz')
parser.add_argument('--theta-right', type=float, default=0.4,
                    help='angle to travel through clockwise, in cycles')
parser.add_argument('--theta-left', type=float, default=0.4,
                    help='angle to travel through counter-clockwise, in cycles')
args = parser.parse_args()

tau = 2.*np.pi
th_right = args.theta_right*tau
th_left = args.theta_left*tau

# Model S can be downloaded from
# http://astro.phys.au.dk/~jcd/solar_models/fgong.l5bi.d.15c
try:
    fgong = io.load_fgong('data/modelS.fgong')
except IOError:
    import urllib2
    response = urllib2.urlopen('http://astro.phys.au.dk/~jcd/solar_models/fgong.l5bi.d.15c')
    with open('data/modelS.fgong','w') as f:
        f.write(response.read())

    fgong = io.load_fgong('data/modelS.fgong')

omega = args.nu*tau # angular frequency in Hz, after JCD's cover picture
    
M,R = fgong['glob'][:2]
r,P,rho,G1,A = fgong['var'][:-1,[0,3,4,9,14]].T
m = M*np.exp(fgong['var'][:-1,1])
cs2 = G1*P/rho

G = 6.672e-8
g = G*m/r**2
Hp = P/rho/g

omega_BV2 = A*g/r
omega_AC2 = cs2/4./Hp**2

t = np.linspace(0., 10000., 10000)/R
x0 = [0.9995*R, tau/4.]

for ell in args.ell:
    k_h = np.sqrt(1.0*ell*(ell+1))/r
    k_r2 = (omega**2-omega_AC2)/cs2 - k_h**2*(1.-omega_BV2/omega**2)
    k_r = np.sqrt(k_r2)

    v_gr = k_r*omega**3*cs2/(omega**4-k_h**2*cs2*omega_BV2)  # dlnr
    v_gh = k_h*omega*cs2*(omega**2-omega_BV2)/(omega**4-k_h**2*cs2*omega_BV2)  # dtheta
    v_gr_k_r = k_r2*omega**3*cs2/(omega**4-k_h**2*cs2*omega_BV2)

    # I = (np.isfinite(v_gr*v_gh))
    # v_gr = v_gr[I]
    # v_gh = v_gh[I]
    # r = r[I]

    # RHS of dlnr/ds, dtheta/ds
    def v(x, t, sign=(-1,-1)):
        ri, thi = x
        # return [sign[0]*ri*np.interp(ri, r[::-1], v_gr[::-1], left=np.nan, right=np.nan),
        return [sign[0]*ri*np.interp(ri, r[::-1], v_gr_k_r[::-1], left=np.nan, right=np.nan) \
                /np.sqrt(np.interp(ri, r[::-1], k_r2[::-1], left=np.nan, right=np.nan)),
                sign[1]*np.interp(ri, r[::-1], v_gh[::-1], left=np.nan, right=np.nan)]

    sol = spint.odeint(v, x0, t, args=((-1,-1),))

    s, th = sol.T
    I = np.isfinite(s*th)
    s = s[I]
    th = th[I]

    # this completes one arc, which we then store
    th = np.hstack((th, 2.*th[-1]-th[::-1]))
    s = np.hstack((s, s[::-1]))

    th_one = np.copy(th)
    s_one = np.copy(s)
    
    # to make the complete arc (with bounces), we add more arcs until
    # we get to the desired angle, then cut everything up to that
    # angle

    # first clockwise
    for i in range(100):
        if np.abs(th[-1]-th[0]) < th_right:
            th = np.hstack((th, th_one-th[0]+th[-1]))
            s = np.hstack((s, s_one))
        else:
            s = s[np.abs(th-th[0]) < th_right]
            th = th[np.abs(th-th[0]) < th_right]
            break
    
    x = s*np.cos(th)/R
    y = s*np.sin(th)/R
    arc, = pl.plot(x, y)

    pl.arrow(x[-2], y[-2], x[-1]-x[-2], y[-1]-y[-2],
             color=arc.get_color(), width=0.002)

    # then counter-clockwise
    th_one = 2.*th_one[0]-th_one
    th = np.copy(th_one)
    s = np.copy(s_one)

    for i in range(100):
        if np.abs(th[-1]-th[0]) < th_right:
            th = np.hstack((th, th_one-th[0]+th[-1]))
            s = np.hstack((s, s_one))
        else:
            s = s[np.abs(th-th[0]) < th_right]
            th = th[np.abs(th-th[0]) < th_right]
            break
    
    x = s*np.cos(th)/R
    y = s*np.sin(th)/R
    arc, = pl.plot(x, y, color=arc.get_color())

    pl.arrow(x[-2], y[-2], x[-1]-x[-2], y[-1]-y[-2],
             color=arc.get_color(), width=0.002)
        
    # then dashed circle for inner turning point
    th = np.linspace(0., 2.*np.pi, 100)
    x = np.min(s)*np.cos(th)/R
    y = np.min(s)*np.sin(th)/R
    pl.plot(x, y, '--', color=arc.get_color())
    
th = np.linspace(0., 2.*np.pi, 100)
s = np.ones(len(th))
x = s*np.cos(th)
y = s*np.sin(th)
    
pl.plot(x, y, 'k-')
pl.subplots_adjust(top=1,bottom=0,left=0,right=1)
pl.axis([-1.1, 1.1, -1.1, 1.1])
pl.axis('off')
pl.show()
