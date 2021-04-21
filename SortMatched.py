import numpy as np
import pandas as pd
from astropy.io import fits
import os
from astropy.table import Table

fname = 'c:/Users/sahal/Desktop/RunawayDetector_4-15a.fits'
data = np.array(fits.open(fname)[1].data)
df = pd.DataFrame(data.byteswap().newbyteorder())

def get_dist(l, b, avg_l, avg_b):
        d = np.sqrt((l-avg_l)**2 + (b-avg_b)**2)
        return d

def getl1(tabl):
        over180 = np.where(tabl>180)[0]
        tabl.iloc[over180] = tabl.iloc[over180]-360
        return tabl

l_projected = df['l'] - (df['vlsrl'] - df['avg_pml'] / 13600000 * np.power(10,df['traceback']))
cond = np.where((df['avg_l'] < 90) | (360 - df['avg_l'] < 90 ))[0]
l_projected.iloc[cond] = df['l1'].iloc[cond] - (df['vlsrl'].iloc[cond]-df['avg_pml'].iloc[cond]) / 13600000 * np.power(10, df['traceback'].iloc[cond])
b_projected = df['b'] - (df['vlsrb']-df['avg_pmb']) / 13600000 * np.power(10, df['traceback'])


dist = np.array(get_dist(l_projected, b_projected, getl1(df['avg_l']), df['avg_b']))

df = df.iloc[np.argsort(dist)]
x = np.unique(df['source_id'], return_index=True)[1]
df = df.iloc[x]

df = df.iloc[np.argsort(df['cluster_label'])]

t = Table.from_pandas(df) 
t.write('c:/users/sahal/desktop/RunawayDetector_4-15a-sorted.fits', overwrite=True)