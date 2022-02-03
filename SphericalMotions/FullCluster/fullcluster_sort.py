### Spherical matching but with all stars in clusters
import argparse

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import LSR, ICRS, SkyCoord
from astropy import units as u

from itertools import combinations, product

from scipy.spatial.distance import squareform

from tqdm import tqdm

kc19_fname = '/Users/aidanmcbride/Documents/Sagitta-Runaways/final6age.fits'
sagitta_fname = '/Users/aidanmcbride/Documents/Sagitta-Runaways/Sagitta_HDBSCAN_named_.fits'
output_fname = '/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/out-yymmdd.fits'

from FullCluster_functions import (maketable, clusterspatialparams, clustermotionparams, getalignment, getregion, get_dist,
                get_relative_velocity, get_tracebacktime, within_distance, bearing)

from FullCluster import sort_table

fname = '/Users/aidanmcbride/Documents/Sagitta-Runaways/proofofconcept.fits'
data = np.array(fits.open(fname)[1].data)
tab = pd.DataFrame(data.byteswap().newbyteorder())


tbt = get_tracebacktime(tab['l'].values, tab['b'].values, tab['vlsrl'].values,  tab['vlsrb'].values, tab['l_0'].values, tab['b_0'].values, tab['vlsrl_0'].values, tab['vlsrb_0'].values)




# sort_table(df.iloc[:int(len(df)/20)])