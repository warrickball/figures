#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl
from tomso import io
from scipy import integrate as spint
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('ell', type=int, nargs='+')
args = parser.parse_args()

tau = 2.*np.pi

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

nu = 3.09e-3  # cyclic frequency in Hz
omega = nu*tau
    
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
    t = t[I]
    x = s*np.cos(th)/R
    y = s*np.sin(th)/R

    # pl.plot(t, r)
    # pl.plot(t, th)
    # pl.show()

    inward_line, = pl.plot(x, y)
    # pl.plot([x[0]], [y[0]], 'o')

    pl.arrow(x[-2], y[-2], x[-1]-x[-2], y[-1]-y[-2],
             color=inward_line.get_color(), width=0.002)

    sol = spint.odeint(v, [s[-1], th[-1]], t, args=((1,-1),))

    s, th = sol.T
    x = s*np.cos(th)/R
    y = s*np.sin(th)/R
    pl.plot(x, y, color=inward_line.get_color())

th = np.linspace(0., 2.*np.pi, 100)
s = np.ones(len(th))
x = s*np.cos(th)
y = s*np.sin(th)
    
pl.plot(x, y, 'k-')
pl.show()

# pl.plot(r, v_gh)
# pl.plot(r, v_gr)
# pl.show()
