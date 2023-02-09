from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('clusters', type=str, nargs='+',
                    help='names of clusters')
parser.add_argument('--list', action='store_true',
                    help="list available clusters")
parser.add_argument('-x', type=str, default='B-V',
                    help="x data: single passband (e.g. V) "
                    "or colour (e.g. B-V) (default=B-V)")
parser.add_argument('-y', type=str, default='V',
                    help="y data: single passband (e.g. V) "
                    "or colour (e.g. B-V) (default=V)")
parser.add_argument('--source', type=str)
parser.add_argument('-z', type=str, default='FeH',
                    help="z (colour) variable (default='FeH')")
args = parser.parse_args()

import numpy as np
import matplotlib.pyplot as pl
from astroquery.vizier import Vizier

# Sarajedini et al. (2007)
# https://ui.adsabs.harvard.edu/abs/2007AJ....133.1658S
# 99% sure V and I are actually F606W and F814W
# catalog = 'J/AJ/133/1658'

# Stetson et al. (2019)
# https://ui.adsabs.harvard.edu/abs/2019MNRAS.485.3042S
catalog = 'J/MNRAS/485/3042'

# Harris (1996)
# https://ui.adsabs.harvard.edu/abs/1996AJ....112.1487H
# already cross-matched in Stetson+2019
# catalog = 'VII/202'

v = Vizier(columns=['**', '+_r'], row_limit=-1)

if args.list:
    cat = v.get_catalogs(catalog+'/table2')[0]
    print(' '.join([k for k in cat['Cluster']]))
    exit(0)

q = v.query_object(args.clusters, catalog=catalog)
vmin = q[0]['__Fe_H_'].min()
vmax = q[0]['__Fe_H_'].max()

for i, row in enumerate(q[0]):
    cluster = q[1][q[1]['Cluster'] == row['Cluster']]
    cluster = cluster[cluster['Sharp'] < 0.3] # & (cluster['e_Bmag'] < 0.2) & (cluster['e_Vmag'] < 0.2)]

    data = np.zeros((2,len(cluster)))
    for i, arg in enumerate([args.x, args.y]):
        band = arg.split('-')
        try:
            data[i] = cluster['%smag' % band[0]] - cluster['%smag' % band[1]]
        except IndexError:
            data[i] = cluster['%smag' % band[0]] - row['__M-m_V']

    pl.scatter(*data, s=5, alpha=0.5, c=[row['__Fe_H_']]*len(cluster), vmin=vmin, vmax=vmax)
        
pl.xlabel(args.x)
pl.ylabel('M_'+args.y)
pl.title('NGC 288 (Stetson et al. 2019)')
pl.gca().invert_yaxis()
pl.show()
