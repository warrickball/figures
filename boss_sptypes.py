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
    
# sptypes = ['%s%i' % (k,i) for k in 'OBAFGKML' for i in range(1,10)]
print(sptypes)

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
        # try:
        #     response = urllib2.urlopen(url % (k+'_+0.0_Dwarf'))
        # except urllib2.HTTPError:
        #     try:
        #         response = urllib2.urlopen(url % k)
        #     except urllib2.HTTPError:
        #         continue

        with open('data/tmp.fits','w') as f:
            f.write(response.read())

        data.append(np.copy(fits.open('data/tmp.fits')[1].data))

    data = np.array(data)
    np.save('data/boss_sptypes.npy', data)

print(data[0].dtype.names)
# pl.plot(data[0]['LogLam'], data[0]['Flux'])
# for i in data:
#     pl.plot(i['LogLam'], i['Flux'])
    
# pl.show()

x = data[0]['LogLam']
y = np.vstack([i['Flux'] for i in data])
y = y/np.max(y, axis=1).reshape((-1,1))
z = pl.cm.jet(y*0.+(x-np.min(x))/(np.max(x)-np.min(x)))
z = 0.1+0.9*z
z[...,-1] = y

kwargs = {'aspect':'auto', 
          'extent':[x[0],x[-1],0,1]}

pl.imshow(np.zeros(z.shape)[:,:,:3], **kwargs)
pl.imshow(z, **kwargs)
pl.grid(False)
pl.show()
