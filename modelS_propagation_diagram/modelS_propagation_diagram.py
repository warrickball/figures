#!/usr/bin/env python

import numpy as np
from tomso import io
from matplotlib import pyplot as pl

# Model S can be downloaded from
# http://astro.phys.au.dk/~jcd/solar_models/fgong.l5bi.d.15c
try:
    fgong = io.load_fgong('modelS.fgong')
except IOError:
    import urllib2
    response = urllib2.urlopen('http://astro.phys.au.dk/~jcd/solar_models/fgong.l5bi.d.15c')
    with open('modelS.fgong','w') as f:
        f.write(response.read())

    fgong = io.load_fgong('modelS.fgong')
        
M,R = fgong['glob'][:2]
r,P,rho,G1,A = fgong['var'][:-1,[0,3,4,9,14]].T
m = M*np.exp(fgong['var'][:-1,1])
cs2 = G1*P/rho

G = 6.672e-8

N = np.sqrt(G*m/r**3*A)/2./np.pi
N[np.isinf(N)] = np.nan
N[np.isnan(N)] = 0.
N[N<=0.] = 1e-99
S = [np.sqrt(1.*l*(l+1)*cs2/r**2)/2./np.pi for l in [1,20]]

pl.plot([-1,2], [5.2e3, 5.2e3], 'k--')
pl.semilogy(r/R, N*1e6)
for Si in S:
    pl.semilogy(r/R, Si*1e6)
    
pl.axis([0., np.max(r/R), 1e1, 1e4])
pl.xlabel('radius (solar units)')
pl.ylabel('cyclic frequency ($\mu$Hz)')
pl.show()
