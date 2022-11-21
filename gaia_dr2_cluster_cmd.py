#!/usr/bin/env python3

# reproduces open cluster CMD (Fig. 2) from Gaia DR2: Observational HRDs
# https://ui.adsabs.harvard.edu/abs/2018A%26A...616A..10G
# https://cdsarc.unistra.fr/viz-bin/cat/J/A+A/616/A10

import numpy as np
import matplotlib.pyplot as pl

# data for open clusters copied by hand from Table 2
ocs = np.genfromtxt('data/gaia_dr2_ocs.dat', names=True, encoding='utf-8', dtype=None, delimiter=',')

try:
    # try to load cached data
    data = np.load('data/gaia_dr2_cluster_cmd.npy', allow_pickle=True)
except:
    # if we fail, download the data again
    from astropy.io.votable import from_table, writeto
    from astropy.table import vstack
    from astroquery.vizier import Vizier
    from astroquery.gaia import Gaia

    # retrieve lists of stars from Vizier
    v = Vizier(columns=['**', '+_r'], row_limit=-1)
    cats = v.get_catalogs('J/A+A/616/A10')

    # tablea1a.dat      71     5378   Stars in nine open clusters within 250pc
    # tablea1b.dat      58    35525   Stars in 37 open clusters beyond 250pc
    # tablea3.dat      160        9   Mean parameters for clusters within 250pc
    # tablea4.dat      120       37   Mean parameters for clusters beyond 250pc

    sources = vstack([cats[0], cats[1]])

    # download Gaia data synchronously, 1900 rows at a time
    # by writing source names to a local VOTable
    # then uploading it and merging against main Gaia DR2 source table
    query = """SELECT *
    FROM tap_upload.table_test AS sources
    INNER JOIN gaiadr2.gaia_source AS gaia
    ON sources.source = gaia.source_id"""

    upload_resource = '/tmp/gaia_dr2_cluster_cmd.xml'

    def get_one(start, end):
        writeto(from_table(sources[start:end], table_id='cluster_cmd'), upload_resource)
        return Gaia.launch_job(query=query, upload_resource=upload_resource,
                               upload_table_name="table_test", verbose=True).get_results()

    data = vstack([get_one(i, i+1900) for i in range(0, len(sources), 1900)]).as_array()

    # cache file locally
    np.save('data/gaia_dr2_cluster_cmd.npy', data, allow_pickle=True)

vmin = ocs['age'].min()
vmax = ocs['age'].max()
for row in ocs[np.argsort(ocs['age'])]:
    cluster = row['Cluster']
    I = data['cluster'] == cluster
    pl.scatter(data[I]['bp_rp'], data[I]['phot_g_mean_mag']-row['DM'], s=3,
               color=pl.cm.jet((row['age']-vmin)/(vmax-vmin)), label=cluster)

pl.xlabel(r"$\mathrm{G}_\mathrm{BP}-\mathrm{G}_\mathrm{RP}$")
pl.ylabel(r"$M_\mathrm{G}$")
pl.gca().invert_yaxis()
pl.colorbar(pl.cm.ScalarMappable(cmap=pl.cm.jet,
                                 norm=pl.Normalize(vmin=vmin, vmax=vmax)),
            label=r"$\log_{10}(\mathrm{age})$")
pl.show()
