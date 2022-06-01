#!/usr/bin/env python3

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('filters', type=str, nargs='+')
parser.add_argument('--legend', type=str, nargs='+', default=[])
parser.add_argument('--interp', type=str, default='none',
                    choices={'none', 'pchip'})
args = parser.parse_args()

import matplotlib.pyplot as pl
from astroquery.svo_fps import SvoFps

if args.interp == 'pchip':
    import numpy as np
    from scipy.interpolate import PchipInterpolator as interp

for name in args.filters:
    data = SvoFps.get_transmission_data(name)
    if args.interp != 'none':
        x = np.linspace(data['Wavelength'].min(),
                        data['Wavelength'].max(),
                        len(data)*10)
        y = interp(data['Wavelength']/10, data['Transmission'])(x/10)
        pl.plot(x, y)
    else:
        pl.plot(data['Wavelength']/10, data['Transmission'])

if len(args.legend) > 0:
    pl.legend(args.legend)
else:
    pl.legend([filter.split('/')[-1] for filter in args.filters])

pl.xlabel('wavelength (nm)')
pl.ylabel('transmission fraction')
pl.ylim(0, 1.05)
pl.show()
