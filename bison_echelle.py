#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl

def get(source, target):
    try:
        data = np.loadtxt(target)
    except IOError:
        try:
            from urllib2 import urlopen
        except ImportError:
            from urllib.request import urlopen

        response = urlopen(source)
        with open(target, 'wb') as f:
            f.write(response.read())

        data = np.loadtxt(target)

    return data

# Davies et al. (2014), doi:http://mnrasl.oxfordjournals.org/content/396/1/L100
# http://bison.ph.bham.ac.uk/downloads/data/davies2014.txt
# Broomhall et al. (2009), doi:10.1007/s11207-015-0810-0
# http://bison.ph.bham.ac.uk/downloads/data/broomhall2009.txt
broomhall = get('http://bison.ph.bham.ac.uk/downloads/data/broomhall2009.txt',
                'data/broomhall2009.txt')
davies = get('http://bison.ph.bham.ac.uk/downloads/data/davies2014.txt',
             'data/davies2014.txt')

Delta_nu = 135.1

data = {}
for row in broomhall:
    data[int(row[0]), int(row[1])] = np.array([row[2], np.mod(row[2], Delta_nu), row[3]])

for row in davies:
    data[int(row[0]), int(row[1])] = np.array([row[2], np.mod(row[2], Delta_nu), row[3]])

styles = 'ovsd'
for ell, style in enumerate(styles):
    nu = np.array([v[0] for k,v in data.items() if k[1]==ell])
    pl.plot(np.mod(nu, Delta_nu), nu, style,
            label=r'$\ell=%i$' % ell)

a = np.array(pl.axis())
pl.axis([0., Delta_nu, a[2], a[3]])
pl.ylabel(r'frequency ($\mu$Hz)')
pl.xlabel(r'frequency module large separation ($\mu$Hz)')
pl.legend(loc=6)
pl.show()

