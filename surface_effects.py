#!/usr/bin/env python

import numpy as np
from tomso import adipls
from matplotlib import pyplot as pl
from scipy.optimize import leastsq
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--correction', type=str, nargs='+',
                    default=[])
args = parser.parse_args()

def get(source, target):
    try:
        data = np.loadtxt(target)
    except IOError:
        try:
            from urllib import urlretrieve
        except ImportError:
            from urllib.request import urlretrieve

        urlretrieve(source, target)
        data = np.loadtxt(target)

    return data

css, eigs = adipls.load_amde('data/modelS.amde')

# Davies et al. (2014), doi:http://mnrasl.oxfordjournals.org/content/396/1/L100
# http://bison.ph.bham.ac.uk/downloads/data/davies2014.txt
# Broomhall et al. (2009), doi:10.1007/s11207-015-0810-0
# http://bison.ph.bham.ac.uk/downloads/data/broomhall2009.txt
broomhall = get('http://bison.ph.bham.ac.uk/downloads/data/broomhall2009.txt',
                'data/broomhall2009.txt')
davies = get('http://bison.ph.bham.ac.uk/downloads/data/davies2014.txt',
             'data/davies2014.txt')

data = np.vstack((davies, broomhall))
data = data[np.unique(1000*data[:,1] + data[:,0],
                      return_index=True)[1]]
data = data[np.argsort(data[:,2])]
n, l, obs, err = data.T

mdl = np.squeeze([css[(css['enn'] == enn) & (css['ell'] == ell)]['nu_Ri']
                  for (enn, ell) in zip(n,l)])*1e3  # mHz -> uHz
E = np.squeeze([css[(css['enn'] == enn) & (css['ell'] == ell)]['E']
                for (enn, ell) in zip(n,l)])*1e3  # mHz -> uHz

styles = 'ovsd'
for ell, style in enumerate(styles):
    pl.plot(obs[l==ell], (obs-mdl)[l==ell], style,
            label=r'$\ell=%i$' % ell)

for corr in args.correction:
    if corr == 'cubic':
        X = (mdl**3/(E/err)).reshape((-1,1))
        y = (obs-mdl)/err
        a3 = np.linalg.lstsq(X, y, rcond=None)[0]
    elif corr == 'both':
        X = np.power(mdl.reshape((-1,1)), [[-1,3]])/(E*err).reshape((-1,1))
        y = (obs-mdl)/err
        a1, a3 = np.linalg.lstsq(X, y, rcond=None)[0]
    elif corr == 'sonoi':
        def func(z):
            return z[0]*numax*(1.-1./(1.+(mdl/numax)**z[1]))            
        
        Q = np.hstack((np.ones(sum(l==0)),
                       E[l>0]/np.interp(mdl[l>0], mdl[l==0], E[l==0])))
        numax = 3100.0
        p = leastsq(lambda z: (func(z) + Q*(mdl-obs))-err, [-0.002, 10])[0]

# a = np.array(pl.axis())
# pl.axis([0., Delta_nu, a[2], a[3]])
pl.xlabel(r'frequency ($\mu$Hz)')
pl.ylabel(r'frequency difference ($\mu$Hz)')
pl.legend(loc=6)
pl.show()

