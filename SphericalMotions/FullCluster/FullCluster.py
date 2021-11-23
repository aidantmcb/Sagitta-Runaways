### Spherical matching but with all stars in clusters
import argparse

import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import LSR, ICRS, SkyCoord
from astropy import units as u

kc19_fname = '/Users/aidanmcbride/Documents/Sagitta-Runaways/final6age.fits'
sagitta_fname = '/Users/aidanmcbride/Documents/Sagitta-Runaways/Sagitta_HDBSCAN_named_.fits'
output_fname = '/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/out-yymmdd.fits'

from FullCluster_functions import (maketable, clusterspatialparams, clustermotionparams, getalignment, getregion, get_dist,
                get_relative_velocity, get_tracebacktime, within_distance)

## Parameters 
yso_cut = 0.9 #select all YSO > .9
parallax_cut = 1 #select all parallax > 1 mas, or distance < 1000 pc


def main():
    sagitta = maketable(sagitta_fname)
    sagitta.rename(columns = {'pms': 'yso'})
    sagitta = sagitta.iloc(np.where((sagitta['yso'] > yso_cut) & (sagitta['parallax'] > parallax_cut))[0]

    #Define l1 column for interactions that may happen across the l = 0 discontinuity
    if 'l1' in sagitta.columns:
        sagitta.drop(columns=['l1']) #reset this column if it already exists
    sagitta['l1'] = sagitta['l']
    over180 = np.where(sagitta['l1']>180)[0]
    sagitta['l1'].iloc[over180] = sagitta['l1'].iloc[over180]-360

    # Extract all labels from clustering, remove the unclustered (-1) sources
    labels = np.sort(np.unique(sagitta['labels']))
    if -1 in labels:
        labels = np.delete(labels, 0)
    print(len(labels))

if __name__ == '__main__':
    main()

def parse_args():
    parser = argparse.Argument_Parser()
        # Main Pipeline Control Options:
    options = parser.add_argument_group("options")
    options.add_argument("tableIn",
                        help="File path and name for table with Gaia source ids (.fits file)",
                        type=str)
    return parser.parse_args()