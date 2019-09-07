#!/usr/bin/env python

import numpy as np
from tomso import fgong
from matplotlib import pyplot as pl

# Model S can be downloaded from
# http://astro.phys.au.dk/~jcd/solar_models/fgong.l5bi.d.15c
try:
    glob, var = fgong.load_fgong('data/modelS.fgong')
except IOError:
    try:
        from urllib2 import urlopen
    except ImportError:
        from urllib.request import urlopen
        
    response = urlopen('http://astro.phys.au.dk/~jcd/solar_models/fgong.l5bi.d.15c')
    with open('data/modelS.fgong','wb') as f:
        f.write(response.read())

    response.close()

    glob, var = fgong.load_fgong('data/modelS.fgong')
        
M, R = glob[:2]
r, P, rho, G1, A = var[:-1,[0,3,4,9,14]].T
m = M*np.exp(var[:-1,1])
cs2 = G1*P/rho

G = 6.672e-8

N2 = G*m/r**3*A
N2[N2<0.] = 1e-99  # prevents NaNs from sqrt(-1)
N = np.sqrt(N2)/2./np.pi
S = [np.sqrt(1.*l*(l+1)*cs2/r**2)/2./np.pi for l in [1,15]]

pl.plot([-1,2], [5.2e3, 5.2e3], 'k--')
line, = pl.semilogy(r/R, N*1e6)
pl.fill_between(r/R, 1e-99, N*1e6, facecolor=line.get_color(), alpha=0.5)
for Si in S:
    line, = pl.semilogy(r/R, Si*1e6)
    pl.fill_between(r/R, Si*1e6, 1e99, facecolor=line.get_color(), alpha=0.5)

pl.text(0.35, 10.**3., '$\ell=1$ pressure modes',
            horizontalalignment='center',
            verticalalignment='bottom')
pl.text(0.8, 10.**3.33, '$\ell=15$ pressure modes',
            horizontalalignment='center',
            verticalalignment='center')
pl.text(0.35, 10.**1.75, 'all gravity modes',
            horizontalalignment='center',
            verticalalignment='center')
pl.text(0.99, 5.5e3, 'acoustic cut-off frequency',
        horizontalalignment='right',
        verticalalignment='bottom')
        
    
pl.axis([0., np.max(r/R), 1e1, 1e4])
pl.xlabel('radius (solar units)')
pl.ylabel('cyclic frequency ($\mu$Hz)')
pl.show()
