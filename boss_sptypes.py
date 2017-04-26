#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot as pl

# Kesseli et al. (2017), arXiv:1702.06957
sptypes = ['O5','O7','O8','O9'] \
          + ['%s%i' % (k,i) for k in 'B' for i in range(1,10)] \
          + ['A0','A1','A2'] + ['A%i_+0.0_Dwarf' % i for i in range(3,10)] \
          + ['%s%i_+0.0_Dwarf' % (k,i) for k in 'FGKM' for i in range(1,10)]
sptypes[-1] = 'M9'
sptypes += ['L0','L1','L2','L3','L6']
sptypes.remove('B7')
for i in [4,5,8]:
    sptypes.remove('A%i_+0.0_Dwarf' % i)

for i in [6,8,9]:
    sptypes.remove('K%i_+0.0_Dwarf' % i)
    
try:
    data = np.load('data/boss_sptypes.npy')
except IOError:
    from astropy.io import fits
    import urllib2
    data = []
    url = 'https://raw.github.com/BU-hammerTeam/PyHammer/master/resources/templates/%s.fits'
    
    for k in sptypes:
        print('Downloading spectrum for SpType %s...' % k)
        response = urllib2.urlopen(url % k)

        with open('data/tmp.fits','w') as f:
            f.write(response.read())

        data.append(np.copy(fits.open('data/tmp.fits')[1].data))

    data = np.array(data)
    np.save('data/boss_sptypes.npy', data)

xx = 10.**data[0]['LogLam']/10.
x = np.linspace(xx[0], xx[-1], len(xx))
y = np.vstack([np.interp(x, xx, i['Flux']) for i in data])
y = y/np.max(y, axis=1).reshape((-1,1))
z = pl.cm.jet(y*0.+(x-np.min(x))/(np.max(x)-np.min(x)))
z = 0.1+0.9*z
z[...,-1] = y

kwargs = {'aspect':'auto', 
          'extent':[x[0],x[-1],0,1]}

# show rainbow of spectra
pl.imshow(np.zeros(z.shape)[:,:,:3], **kwargs)
pl.imshow(z, **kwargs)
pl.xlabel('wavelength (nm)')
pl.yticks([])
pl.grid(False)

# annotate spectral lines
# line data from http://www.star.ucl.ac.uk/~msw/lines.html
line_data = [['Ha', 656.461],
             ['Hb', 486.269],
             ['NaI', (589.1583+589.7558)/2.],
             ['MgI', 517.4125],
             ['CaII', 854.444]]

A = pl.axis()
for (label, wavelength) in line_data:
    label = label.replace('Ha',r'H$\alpha$')
    label = label.replace('Hb',r'H$\beta$')
    line, = pl.plot([wavelength,wavelength],[0,1], 'w-', alpha=0.25, lw=1.5);
    line.set_dashes([4,6])
    pl.annotate(s=label, xy=((wavelength-A[0])/(A[1]-A[0]), 1.02),
                xycoords='axes fraction', color='k', ha='center')

# annotate spectral types
N = len(sptypes)
i_prev = 0
yticks = [1.0]
for letter in ['O','B','A','F','G','K','M','L']:
    i_next = [i for i in range(N) if sptypes[i].startswith(letter)][-1]+1
    yticks.append(1.0-1.0*i_next/N)
    y_text = 1.0-1.0*i_prev/N-0.5*(i_next-i_prev)/N
    pl.annotate(s=letter,xy=(-0.04, y_text), xycoords='axes fraction',
                va='center')
    i_prev = i_next

pl.yticks(yticks, ['' for i in range(len(yticks))])
pl.annotate(s='spectral type', xy=(-0.08, 0.5), xycoords='axes fraction',
            va='center', rotation=90.0)
pl.show()
