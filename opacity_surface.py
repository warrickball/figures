#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-X', type=str, default='0.02',
                    help="string for Z (default='0.02')")
parser.add_argument('-Z', type=str, default='0.7',
                    help="string for X (default='0.7')")
args = parser.parse_args()

key = 'z%s_x%s' % (args.Z, args.X)

import numpy as np
from matplotlib import pyplot as pl
import os

with open('%s/data/kap_data/gs98_%s.data' % (os.environ['MESA_DIR'], key),'r') as f:
    lines = f.readlines()

opac1 = np.loadtxt(lines[7:])[:,1:].T
logT1 = np.loadtxt(lines[7:])[:,0]
logR1 = np.loadtxt(lines[5:6])

logrho1 = logR1[:,np.newaxis]+3.*logT1-18.
logT1 = 0.*logR1[:,np.newaxis]+logT1
i1 = np.where(logT1[0]>=3.85)[0][0]-1

opac1 = opac1[:,i1:]
logrho1 = logrho1[:,i1:]
logT1 = logT1[:,i1:]

with open('%s/data/kap_data/lowT_fa05_gs98_%s.data' % (os.environ['MESA_DIR'], key),'r') as f:
    lines = f.readlines()

opac2 = np.loadtxt(lines[7:])[:,1:].T
logT2 = np.loadtxt(lines[7:])[:,0]
logR2 = np.loadtxt(lines[5:6])

logrho2 = logR2[:,np.newaxis]+3.*logT2-18.
logT2 = 0.*logR2[:,np.newaxis]+logT2
i2 = np.where(logT2[0]<=3.85)[0][-1]

opac2 = opac2[:,:i2]
logrho2 = logrho2[:,:i2]
logT2 = logT2[:,:i2]

levels = np.linspace(min(np.min(opac1),np.min(opac2)),
                     max(np.max(opac1),np.max(opac2)),40)

pl.contourf(logrho1, logT1, opac1, levels)
pl.contourf(logrho2, logT2, opac2, levels)
c = pl.contour(logrho1, logT1, opac1, range(0,6),colors='k')
pl.clabel(c, inline=1, fontsize=8,fmt='%i')
c = pl.contour(logrho1, logT1, opac1, range(-5,0),colors='k')
c = pl.contour(logrho2, logT2, opac2, range(-5,0),colors='k')
pl.clabel(c, inline=1, fontsize=8,fmt='%i')
c = pl.contour(logrho2, logT2, opac2, range(0,6),colors='k')
pl.xlabel(r'$\log_{10}(\rho/(\mathrm{g}/\mathrm{cm}^3))$')
pl.ylabel(r'$\log_{10}(T/\mathrm{K})$')
pl.show()
