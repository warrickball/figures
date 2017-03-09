#!/usr/bin/env python

from matplotlib import pyplot as pl
import numpy as np
from scipy.interpolate import interp1d

hip_dtype = [('Catalog', '|S1'), ('HIP', int), ('Proxy', '|S1'),
             ('RAhms', '|S11'), ('DEdms', '|S11'), ('Vmag', float),
             ('VarFlag', int), ('r_Vmag', '|S1'), ('RAdeg', float),
             ('DEdeg', float), ('AstroRef', '|S1'), ('Plx', float),
             ('pmRA', float), ('pmDE', float), ('e_RAdeg', float),
             ('e_DEdeg', float), ('e_Plx', float), ('e_pmRA', float),
             ('e_pmDE', float), ('DE:RA', float), ('Plx:RA', float),
             ('Plx:DE', float), ('pmRA:RA', float),
             ('pmRA:DE', float), ('pmRA:Plx', float),
             ('pmDE:RA', float), ('pmDE:DE', float),
             ('pmDE:Plx', float), ('mDE:pmRA', float), ('F1', int),
             ('F2', float), ('---', int),
             ('BTmag', float), ('e_BTmag', float), ('VTmag', float),
             ('e_VTmag', float), ('m_BTmag', '|S1'), ('B-V', float),
             ('e_B-V', float), ('r_B-V', '|S1'), ('V-I', float),
             ('e_V-I', float), ('r_V-I', '|S1'), ('CombMag', '|S1'),
             ('Hpmag', float), ('e_Hpmag', float), ('Hpscat', float),
             ('o_Hpmag', int), ('m_Hpmag', '|S1'), ('Hpmax', float),
             ('HPmin', float), ('Period', float), ('HvarType', '|S1'),
             ('moreVar', '|S1'), ('orePhoto', '|S1'),
             ('CCDM', '|S10'), ('n_CCDM', '|S1'), ('Nsys', int),
             ('Ncomp', int), ('MultFlag', '|S1'), ('Source', '|S1'),
             ('Qual', '|S1'), ('m_HIP', '|S2'), ('theta', int), ('rho', float),
             ('e_rho', float), ('dHp', float), ('e_dHp', float),
             ('Survey', '|S1'), ('Chart', '|S1'), ('Notes', '|S1'),
             ('HD', int), ('BD', '|S10'), ('CoD', '|S10'),
             ('CPD', '|S10'), ('(V-I)red', float), ('SpType', '|S12'),
             ('r_SpType', '|S1')]
            
try:
    data = np.load('data/hip_main.npy')
except IOError:
    import urllib2
    response = urllib2.urlopen('ftp://cdsarc.u-strasbg.fr/pub/cats/I/239/hip_main.dat.gz')
    with open('data/hip_main.dat.gz', 'w') as f:
        f.write(response.read())

    np.save('data/hip_main.npy', np.genfromtxt('data/hip_main.dat.gz',
                                               delimiter='|',
                                               dtype=hip_dtype))
    data = np.load('data/hip_main.npy')


# print(data.dtype.names)
for k in ['e_Hpmag', 'e_Plx']:
    data = data[~np.isnan(data[k])]

data = data[data['e_Hpmag']<0.1]
data = data[data['e_Plx']/data['Plx']<0.05]
data = data[data['Plx']>0]

Hp = data['Hpmag'] + 5.0 + 5.0*np.log10(1e-3*data['Plx'])
BV = data['BV']

# original annotations
MS_BV, MS_Hp = np.array([[-0.18,-2.46], [0.0,0.93], [0.2,2.45],
                         [0.4,3.3], [0.63,5.], 
                         [0.9,6.4], [1.2,7.8],
                         [1.37,8.48], [1.57,11.75]]).T
MS_BVi = np.linspace(np.min(MS_BV), np.max(MS_BV), 100)
pl.plot(MS_BVi, interp1d(MS_BV, MS_Hp, kind='cubic')(MS_BVi), 'g-', lw=10,
        alpha=0.6, solid_capstyle='round')  # MS
pl.fill_between([0.25,0.75,0.85],[0.,2.8,2.0], [0.,0.,2.], alpha=0.25,
                lw=0)  # Hertzsprung gap
pl.plot([0.85,1.6], [3.6,-1.05], 'r-', lw=25, alpha=0.2,
        solid_capstyle='round')  # RGB
pl.plot([0.92,1.07],[0.9,1.],'r-', lw=25, alpha=0.2,
        solid_capstyle='round')  # RC
pl.plot([-0.5,1.0], [8.6,17.34], 'b-', lw=25, alpha=0.1)  # WD
pl.plot(0.653, 4.86, 'y*', ms=15)  # Sun

# pl.hexbin(BV, Hp)
pl.plot(BV, Hp, 'ko', ms=4, alpha=0.1)
pl.axis([-0.5, 2., 15., -5.])
pl.xlabel(r'$B-V$')
pl.ylabel(r'$M_\mathrm{Hipparcos}$')
pl.show()

