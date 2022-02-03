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

## Parameters 
yso_cut = 0.9 #select all YSO > .9
parallax_cut = 1 #select all parallax > 1 mas, or distance < 1000 pc

def thing(row1, row2):
    l, b, vlsrl, vlsrb = (0,1,2,3)

    starl = row2[0]
    
    if 360 - starl < 90:
        l = 4
        row2[0] = starl - 360
    elif starl < 90:
        l = 4

    lon3 = (row1[l] - (row1[vlsrl] - row2[vlsrl])/1e3/3600)
    lat3 = (row1[b] - (row1[vlsrb] - row2[vlsrb])/1e3/3600)

    bearing_to = bearing(row1[b], row1[l], row2[b], row2[l])
    bearing_of = bearing(row1[b], row1[l], lat3, lon3)
    alignment = np.cos(bearing_of - bearing_to)

    return alignment


def main():
    sagitta = maketable(sagitta_fname)
    sagitta.rename(columns = {'pms': 'yso'})

    # for column in sagitta.columns:
    #     if sagitta[column].dtype == 'object':
    #         first = sagitta[column].iloc[0]
    #         first = str(first)
    #         if any(c.isalpha for c in first):
    #             sagitta[column] = sagitta[column].astype('str')
    #         elif '.' in first:
    #             sagitta[column] = sagitta['column'].astype(float)
    #         else:
    #             sagitta[column] = sagitta['column'].astype(int)



    sagitta = sagitta.iloc[np.where((sagitta['yso'] > yso_cut) & (sagitta['parallax'] > parallax_cut))[0]]

    #Define l1 column for interactions that may happen across the l = 0 discontinuity
    if ('l1' in sagitta.columns):
        sagitta.drop(columns=['l1']) #reset this column if it already exists
    sagitta['l1'] = sagitta['l']
    over180 = np.where(sagitta['l1']>180)[0]
    sagitta['l1'].iloc[over180] = sagitta['l1'].iloc[over180]-360

    # Extract all labels from clustering, remove the unclustered (-1) sources
    labels = np.sort(np.unique(sagitta['labels']))
    if -1 in labels:
        labels = np.delete(labels, 0)
    print('Labeled Clusters:', len(labels))


    # np.random.seed(1)
    # downsample_by = 50
    # sagitta = sagitta.iloc[np.random.choice(len(sagitta), int(len(sagitta)/downsample_by))]
    # print(len(sagitta))

    d1, d2 = (sagitta.iloc[np.where(sagitta['labels']!=-1)[0]], sagitta.iloc[np.where(sagitta['labels']==-1)[0]])

    # alignments = np.loadtxt('/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/alignments_full.gz')
    # # print(len(alignments))
    # alignments = alignments.reshape((len(d1), len(d2))) #make to a matrix with rows and cols representing indices

    # points1 = np.where(1 - np.abs(alignments) < 0.1)[0]
    # points2 = np.where(1-np.abs(alignments) < 0.1)[1]
    # np.savetxt('/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/ind1.gz', points1)
    # np.savetxt('/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/ind2.gz', points2)
    # print(alignments.shape)
    # print(len(d1) * len(d2))

    points1 = np.loadtxt('/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/ind1.txt')
    points2 = np.loadtxt('/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/ind2.txt')
    
    np.savetxt('/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/ind1_.txt', points1.astype(int))
    np.savetxt('/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/ind2_.txt', points2.astype(int))

    return

    stars1 = sagitta.iloc[points1]
    stars2 = sagitta.iloc[points2]
    
    l = 'l'

    lon3=(stars1[l]-(stars1['vlsrl']-stars2[vlsrl])/1e3/3600)#*np.pi/180
    # lon3=(tab['l']-(tab['vlsrl']-avg_pml)/1e3/3600)#*np.pi/180

    lat3=(stars1['b']-(stars1['vlsrb']-stars2['vlsrb'])/1e3/3600)#*np.pi/180


    bearing_to = bearing(stars1['b'], stars1[l], stars2['b'], stars2[l])
    bearing_of = bearing(stars1['b'], stars1[l], lat3, lon3)
    alignment = np.cos(bearing_of - bearing_to)

    tbt = get_tracebacktime(stars1[l], stars1['b'], stars1['vlsrl'],  stars1['vlsrb'], stars2[l], stars2['b'], stars2['vlsrl'], stars2['vlsrb']) 

    criteria = lambda x, y: (1-x) * 10 + y

    # best_both = minimize(criteria, df['alignment'], df['traceback']/df['avg_age'])


    age = star2['age']
    # df = df.iloc[np.argsort(criteria(alignment, tbt / age))]

    sortthing = np.argsort(criteria(alignment, tbt / age))
    stars1 = stars1.iloc[sortthing]
    stars2 = stars2.iloc[sortthing]


    x = np.unique(stars1["source_id"], return_index=True)[1]
    bestmatch = np.zeros(len(stars1), dtype="bool")
    bestmatch[x] = True
    print('unlikely, it worked')
    # df = df.iloc[x]
    df["best"] = bestmatch
    df = df.iloc[np.argsort(df["cluster_label"])]
    print("Final length: " + str(len(df)))

    t = Table.from_pandas(df)
    t.write(
        "/Users/aidanmcbride/Documents/Sagitta-Runaways/Outputs/RunawayDetector_11-5-21-AlignTimeSorted.fits",
        overwrite=True,
    )

main()