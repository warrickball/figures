#!/usr/bin/env python

"""Illustrates (approximate) ray propagation for p-mode in Model S.
Currently only seems to work up to l=49 (and fails on a few
intermediate values too).  Still need to
- check that the inner turning point is correct,
- start integration from upper reflecting boundary, and
- add options for how many times to "bounce" at the surface.

"""

# formulae from http://soi.stanford.edu/papers/dissertations/giles/thesis/PDF/chapter02.pdf

import numpy as np
from matplotlib import pyplot as pl
from tomso import fgong
from scipy import integrate as spint
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-l', '--ell', type=int, nargs='+', default=[2, 20, 25, 75],
                    help="angular degree(s) of the desired rays "
                    "(default 2, 20, 25, 75)")
parser.add_argument('-f', '--freq', type=float, default=3.0,
                    help="cyclic frequency in mHz (default=3.0, as in cover of "
                    "JCD's notes")
parser.add_argument('-o', '--output', type=str,
                    help="save figure to file instead of plotting")
parser.add_argument('--theta-right', type=float, default=0.4,
                    help="angle to travel through clockwise, in cycles "
                    "(default=0.4)")
parser.add_argument('--theta-left', type=float, default=0.4,
                    help="angle to travel through counter-clockwise, in cycles "
                    "(default=0.4)")
parser.add_argument('--figsize', type=float, nargs=2,
                    help="figure size, passed to rcParams['figure.figsize']")
parser.add_argument('--padding', type=float, default=0.01,
                    help="fractional padding between edge and circle "
                    "(default=0.01)")
args = parser.parse_args()

if args.figsize:
    pl.rcParams['figure.figsize'] = args.figsize

TAU = 2.*np.pi
th_right = args.theta_right*TAU
th_left = args.theta_left*TAU

# Model S can be downloaded from
# http://astro.phys.au.dk/~jcd/solar_models/fgong.l5bi.d.15c
try:
    S = fgong.load_fgong('data/modelS.fgong', G=6.67232e-8, return_object=True)
except IOError:
    try:
        from urllib2 import urlopen
    except ImportError:
        from urllib.request import urlopen
        
    response = urlopen('http://astro.phys.au.dk/~jcd/solar_models/fgong.l5bi.d.15c')
    with open('data/modelS.fgong','wb') as f:
        f.write(response.read())

    response.close()

    S = fgong.load_fgong('data/modelS.fgong', G=6.67232e-8, return_object=True)

omega = args.freq*TAU*1e-3

omega_AC2 = S.cs2/4./S.Hp**2

t = np.linspace(0., 10000., 100000)/S.R
x0 = [0.9995*S.R, TAU/4.]

for ell in args.ell:
    k_h = np.sqrt(1.0*ell*(ell+1))/S.r
    k_r2 = (omega**2-omega_AC2)/S.cs2 - k_h**2*(1.-S.N2/omega**2)
    k_r = np.sqrt(k_r2)

    v_gr = k_r*omega**3*S.cs2/(omega**4-k_h**2*S.cs2*S.N2)  # dlnr
    v_gh = k_h*omega*S.cs2*(omega**2-S.N2)/(omega**4-k_h**2*S.cs2*S.N2)  # dtheta
    v_gr_k_r = k_r2*omega**3*S.cs2/(omega**4-k_h**2*S.cs2*S.N2)

    # I = (np.isfinite(v_gr*v_gh))
    # v_gr = v_gr[I]
    # v_gh = v_gh[I]
    # r = r[I]

    # RHS of dlnr/ds, dtheta/ds
    def v(x, t, sign=(-1,-1)):
        ri, thi = x
        # return [sign[0]*ri*np.interp(ri, r[::-1], v_gr[::-1], left=np.nan, right=np.nan),
        return [sign[0]*ri*np.interp(ri, S.r[::-1], v_gr_k_r[::-1], left=np.nan, right=np.nan) \
                /np.sqrt(np.interp(ri, S.r[::-1], k_r2[::-1], left=np.nan, right=np.nan)),
                sign[1]*np.interp(ri, S.r[::-1], v_gh[::-1], left=np.nan, right=np.nan)]

    sol = spint.odeint(v, x0, t, args=((-1,-1),))

    s, th = sol.T
    I = np.isfinite(s*th)
    s = s[I]
    th = th[I]

    # truncates solutions that accidentally converge on dtheta/dr=0
    I = np.where(np.diff(th)/np.diff(s)<0.)[0]
    if len(I)>0:
        s = s[:I[0]]
        th = th[:I[0]]
    
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
    
    x = s*np.cos(th)/S.R
    y = s*np.sin(th)/S.R
    arc, = pl.plot(x, y)

    pl.arrow(x[-2], y[-2], x[-1]-x[-2], y[-1]-y[-2],
             color=arc.get_color(), width=0.01)

    # then counter-clockwise
    th_one = 2.*th_one[0]-th_one
    th = np.copy(th_one)
    s = np.copy(s_one)

    for i in range(100):
        if np.abs(th[-1]-th[0]) < th_left:
            th = np.hstack((th, th_one-th[0]+th[-1]))
            s = np.hstack((s, s_one))
        else:
            s = s[np.abs(th-th[0]) < th_left]
            th = th[np.abs(th-th[0]) < th_left]
            break
    
    x = s*np.cos(th)/S.R
    y = s*np.sin(th)/S.R
    arc, = pl.plot(x, y, color=arc.get_color())

    pl.arrow(x[-2], y[-2], x[-1]-x[-2], y[-1]-y[-2],
             color=arc.get_color(), width=0.01)
        
    # then dashed circle for inner turning point
    th = np.linspace(0., TAU, 100)
    s_t = np.interp(0., S.cs2-omega**2/ell/(ell+1)*S.r**2, S.r)
    x = s_t*np.cos(th)/S.R
    y = s_t*np.sin(th)/S.R
    pl.plot(x, y, '--', color=arc.get_color())
    
th = np.linspace(0., TAU, 100)
s = np.ones(len(th))
x = s*np.cos(th)
y = s*np.sin(th)
    
pl.plot(x, y, 'k-')
b = args.padding
pl.subplots_adjust(top=1-b, bottom=b, left=b, right=1-b)
pl.axis('off')

if args.output:
    pl.savefig(args.output)
else:
    pl.show()
