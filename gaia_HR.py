#!/usr/bin/env python3

# intended to mimic Gaia DR2 HR diagrams (A&A 616, A10, 2018),
# specifically, Fig. 6c because I can't be bothered with this extinction business
# https://www.aanda.org/articles/aa/abs/2018/08/aa32843-18/aa32843-18.html

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.colors import PowerNorm

try:
    r = np.load('data/gaia_HR.npy')
except IOError:
    from astroquery.gaia import Gaia
    job = Gaia.launch_job_async("""
SELECT
phot_g_mean_mag AS g,
parallax,
e_bp_min_rp_val AS e_bp_rp,
phot_bp_mean_mag AS bp,
phot_rp_mean_mag AS rp
FROM gaiadr2.gaia_source
WHERE parallax_over_error > 10
AND parallax > 5
AND phot_g_mean_flux_over_error > 50
AND phot_rp_mean_flux_over_error > 20
AND phot_bp_mean_flux_over_error > 20
AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2)
AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2)
AND visibility_periods_used > 8
AND astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))""")
    r = job.get_results().as_array().data
    np.save('data/gaia_HR.npy', r)

G = r['g']+5.*np.log10(r['parallax'])-10

pl.plot(r['bp']-r['rp'], G, 'k.', ms=2, alpha=0.2, rasterized=True, zorder=0);
pl.hist2d(r['bp']-r['rp'], G, cmap='magma', bins=200, cmin=10, norm=PowerNorm(gamma=1/3), zorder=1)
pl.gca().invert_yaxis()
pl.xlabel(r'$B_\mathrm{p}-R_\mathrm{p}$')
pl.ylabel(r'$G$')
pl.show()
