#!/usr/bin/env python3

from argparse import ArgumentParser
import matplotlib.pyplot as pl
from astroquery.svo_fps import SvoFps

parser = ArgumentParser()
parser.add_argument('filters', type=str, nargs='+')
parser.add_argument('--legend', type=str, nargs='+', default=[])
args = parser.parse_args()

for name in args.filters:
    data = SvoFps.get_transmission_data(name)
    pl.plot(data['Wavelength']/10, data['Transmission'])

if len(args.legend) > 0:
    pl.legend(args.legend)
else:
    pl.legend(args.filters)

pl.xlabel('wavelength (nm)')
pl.ylabel('transmission fraction')
pl.ylim(0, 1.05)
pl.show()
