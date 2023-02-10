#!/usr/bin/env python3

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('clusters', type=str, nargs='+',
                    help='names of clusters')
parser.add_argument('--list', action='store_true',
                    help="list available clusters")
parser.add_argument('-z', type=str, default='__Fe_H_',
                    help="z (colour) variable (default='__Fe_H_')")
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as pl
from astroquery.vizier import Vizier
from astropy.io import fits

v = Vizier(columns=['**', '+_r'], row_limit=-1)

if args.list:
    data = v.get_catalogs('J/AJ/133/1658/clusters')[0]
    print(data)

data = v.query_object(args.clusters, catalog='J/AJ/133/1658/acssggc')[0]
clusters = list(np.unique(data['Cluster']))
meta = v.query_object(clusters, catalog='VII/202')[0]

Z = np.argsort(meta[args.z])
vmin = np.min(meta[args.z])
vmax = np.max(meta[args.z])

for (cluster, row) in zip([clusters[i] for i in Z], meta[Z]):
    I = np.where(data['Cluster'] == cluster)[0]
    I = I[data[I]['e_Vmag'] < 0.01]
    I = I[data[I]['e_V-I'] < 0.01]

    if len(args.clusters) > 1:
        z = row[args.z]*np.ones(len(data[I]))
        pl.scatter(data[I]['V-I'], data[I]['Vmag']-row['__m-M_V'], s=3, label=cluster,
                   c=z, vmin=vmin, vmax=vmax, cmap='jet')
    else:
        pl.scatter(data[I]['V-I'], data[I]['Vmag']-row['__m-M_V'], s=3, label=cluster)

pl.gca().invert_yaxis()
pl.xlabel(r"$\mathrm{V}-\mathrm{I}$")
pl.ylabel(r"$M_\mathrm{V}$")

if len(args.clusters) > 1:
    pl.colorbar()
    pl.legend()

pl.show()
