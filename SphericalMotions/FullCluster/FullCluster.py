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

# def stacktables(tab, tab1):
    # cols = ['source_id', 'l', 'l1', 'b', 'parallax', 'vlsrl', 'vlsrb', 'labels', 'age']
    # x = np.arange(len(tab))
    # y = np.arange(len(tab1))
    # xx, yy = np.meshgrid(x, y)
    
    # xx = np.ravel(xx[1:]) #wait a second this won't work anymore, oh no
    # yy = np.ravel(yy[1:])

    # attach_tab = tab1[cols].iloc[yy]
    # base_tab = tab.iloc[xx]
    # # attach_tab = tab.iloc[yy]
    # attach_tab.rename(columns = lambda name : name + '_0', inplace = True)
    # # print(attach_tab)
    # return pd.concat([base_tab, attach_tab.set_index(base_tab.index)], axis = 1)

# def aligns(tab, tolerance = 0.1):
#     retain = stacktables(tab.iloc[:1], tab.iloc[:1])
#     retain['alignment'] = np.zeros(len(retain), dtype = float)
#     retain = pd.DataFrame(columns = retain.columns)
#     step = 1
#     for i in np.arange(0, len(tab), step):
#         print(i, 'out of', len(tab))
#         stacked = stacktables(tab.iloc[i:i+step], tab)

#         l = 'l' #FIX
#         l_0 ='l_0' 
        
#         lon3 = (stacked[l] - (stacked['vlsrl'] - stacked['vlsrl_0'])/1e3/3600)
#         lat3 = (stacked['b'] - (stacked['vlsrb'] - stacked['vlsrb_0'])/1e3/3600)

#         bearing_to = bearing(stacked['b'], stacked[l], stacked['b'], stacked[l_0])
#         bearing_of = bearing(stacked['b'], stacked[l], lat3, lon3)
#         alignment = np.cos(bearing_of - bearing_to)

#         stacked['alignment'] = alignment
#         save = stacked.iloc[np.where(1-np.abs(alignment) < tolerance)[0]]
#         retain = pd.concat([retain, save], axis = 0)
#     return retain
#     #use np.meshgrid - x coordinates will be 


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

def sort_table(tab):
        l_0_col = tab['l_0']
        condition1 = 360 - l_0_col < 90
        condition2 = l_0_col < 90
        l = np.zeros(len(tab), dtype = 'object')
        l[:] = 'l'
        l[condition1] = 'l1'
        l[condition2] = 'l1'

        l_0 = np.zeros(len(tab), dtype = 'object')
        l_0[:] = 'l_0'
        l_0[condition1] = 'l1_0'
        
        tbt = get_tracebacktime(tab[l].values, tab['b'].values, tab['vlsrl'].values,  tab['vlsrb'].values, tab[l_0].values, tab['b_0'].values, tab['vlsrl_0'].values, tab['vlsrb_0'].values) 
        # tab['traceback'] = np.log10(tbt)

        # tab = tab.iloc[np.where(tbt <  np.power(10, tab['age_0'].values))[0]]
        # print(len(tab))

        # def criteria(x, y):
        #     return (1 - x) * 10 + y

        # a = tab['alignment']
        # traceback = tab['traceback']
        # avg_age = tab['age_0']
        # # best_both = minimize(criteria, df['alignment'], df['traceback']/df['avg_age'])
        # tab = tab.iloc[np.argsort(criteria(a, traceback / avg_age))]
        # x = np.unique(tab["source_id"], return_index=True)[1]
        # tab = tab.iloc[x]

        tab = tab.iloc[np.argsort(1-np.abs(tab['alignment']))]
        tab = tab.iloc[np.unique(tab['source_id'], return_index  = True)[1]]
        

        return tab
    

def sort_data(d1, d2):
        
        l_0_col = d2[:, 0]
        condition1 = 360 - l_0_col < 90
        condition2 = l_0_col < 90
        l = np.zeros(len(l_0_col), dtype = int)
        l[condition1] = 4
        l[condition2] = 4 #deals with the l1 stuff

        l_0 = np.zeros(len(l_0_col), dtype = int)
        l_0[condition1] = 4
        
        tbt = get_tracebacktime(d1[:,l], d1[:,1], d1[:,2],  d1[:,3], d2[:,l_0], d2[:,1], d2[:,2], d2[:,3]) 
        # tab['traceback'] = np.log10(tbt)

        # tab = tab.iloc[np.where(tbt <  np.power(10, tab['age_0'].values))[0]]
        # print(len(tab))

        # def criteria(x, y):
        #     return (1 - x) * 10 + y

        # a = tab['alignment']
        # traceback = tab['traceback']
        # avg_age = tab['age_0']
        # # best_both = minimize(criteria, df['alignment'], df['traceback']/df['avg_age'])
        # tab = tab.iloc[np.argsort(criteria(a, traceback / avg_age))]
        # x = np.unique(tab["source_id"], return_index=True)[1]
        # tab = tab.iloc[x]

        # tab = tab.iloc[np.argsort(1-np.abs(tab['alignment']))]
        # tab = tab.iloc[np.unique(tab['source_id'], return_index  = True)[1]]
        

        return tbt
    

def main():
    sagitta = maketable(sagitta_fname)
    sagitta.rename(columns = {'pms': 'yso'})

    for column in sagitta.columns:
        if sagitta[column].dtype == 'object':
            first = sagitta[column].iloc[0]
            first = str(first)
            if any(c.isalpha for c in first):
                sagitta[column] = sagitta[column].astype('str')
            elif '.' in first:
                sagitta[column] = sagitta['column'].astype(float)
            else:
                sagitta[column] = sagitta['column'].astype(int)



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

    
    ##downsample
    # downsample_by = 50

    # np.random.seed(1)
    # sagitta = sagitta.iloc[np.random.choice(len(sagitta), int(len(sagitta)/downsample_by))]

    labeled_stars = sagitta.iloc[np.where(sagitta['labels']!=-1)[0]]
    unlabeled_stars = sagitta.iloc[np.where(sagitta['labels']==-1)[0]]

    print('labeled stars', len(labeled_stars))
    print('unlabeled stars', len(unlabeled_stars))
    print('length', len(sagitta))


    #ones using 'l' column
    d1 = unlabeled_stars[['l', 'b', 'vlsrl', 'vlsrb', 'l1']].values
    d2 = labeled_stars[['l', 'b', 'vlsrl', 'vlsrb', 'l1', 'age']].values #could use THESE to sort after doing pair_first, second
    alignments_array = np.zeros(len(d1) * len(d2))
    counter = 0
    for pair in tqdm(product(d1, d2), total = len(d1) * len(d2)):
        a = thing(*pair)
        alignments_array[counter] = a
        counter += 1

    np.savetxt('/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/alignments_full.gz', alignments_array)

#     alignments = alignments_array.reshape((len(d1), len(d2))) #make to a matrix with rows and cols representing indices
    
#     pair_first, pair_second = (np.where(1-np.abs(alignments) < 0.1)[0], np.where(1-np.abs(alignments) < 0.1)[1])
#     # print(len(pair_first), len(pair_second))

#     left = unlabeled_stars.iloc[pair_first]
#     right = labeled_stars[['source_id', 'l', 'b', 'vlsrl', 'vlsrb', 'parallax', 'age', 'yso', 'labels', 'f6_name', 'l1']].iloc[pair_second]
#     right.rename(columns = lambda name : name + '_0', inplace = True)        
#     output_frame = pd.concat([left, right.set_index(left.index)], axis = 1)
#     output_frame = output_frame.reset_index(drop = True)
#     output_frame['alignment'] = alignments[pair_first, pair_second]

#     # output_frame = output_frame.iloc

#     for column in output_frame.columns:
#         if output_frame[column].dtype == 'object':
#             first = output_frame[column].iloc[0]
#             first = str(first)
#             if any(c.isalpha for c in first):
#                 output_frame[column] = output_frame[column].astype('str')
#             elif '.' in first:
#                 output_frame[column] = output_frame['column'].astype(float)
#             else:
#                 output_frame[column] = output_frame['column'].astype(int)

#     output_length = len(output_frame)
#     print('Length Before Sort', len(output_frame))

#     # output_frame = output_frame.iloc[sort_data(d1[pair_first,:], d2[pair_second, :])]
#     # output_frame = sort_table(output_frame)

#     print('Length After Sort', len(output_frame))

# # save all alignments
# # threshold .1
# # run traceback times for that subset, and compare to relative ages
# # THEN sort - sort criteria possibly binned stars per label 

#     output_table = Table.from_pandas(output_frame)
#     # output_table.write('/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/FullCluster/data_out/proofofconcept_' + str(downsample_by) + '.fits', overwrite=True)

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

"""
1 2 3 4 
1 2 3 4

1 2 3 4
4 1 2 3

1234
3412

for loop
df.iloc[index]

itertools
"""

# tab with length 4
# index = [1,2,3,4]
# for i:
# funct(tab[index], tab[index->i])

# a = [1, 3 ,5 ]
# b = [2 ,4, 6]
# meshgrid(a, b)
# > [[1, 1 , 1], [2 ,4 ,6]] [3 ,3 ,3 ] [2, 4 ,6]