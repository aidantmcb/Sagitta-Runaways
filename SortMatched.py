import numpy as np
import pandas as pd
from astropy.io import fits
import os
from astropy.table import Table

fname = 'c:/Users/sahal/Desktop/RunawayDetector_3-19a.fits'
data = np.array(fits.open(fname)[1].data)
df = pd.DataFrame(data.byteswap().newbyteorder())

def get_dist(l, b, avg_l, avg_b):
        d = np.sqrt((l-avg_l)**2 + (b-avg_b)**2)
        return d

l_projected = df['l'] - (df['vlsrl']-df['avg_pml']) / 13600000 * np.power(10, df['traceback'])
b_projected = df['b'] - (df['vlsrb']-df['avg_pmb']) / 13600000 * np.power(10, df['traceback'])

dist = np.array(get_dist(l_projected, b_projected, df['avg_l'], df['avg_b']))

x = np.unique(df['source_id'].iloc[np.argsort(dist)], return_index=True)[1]
df = df.iloc[x]

df = df.iloc[np.argsort(df['cluster_label'])]

t = Table.from_pandas(df)
t.write('c:/users/sahal/desktop/RunawayDetector_3-19a-sorted.fits', overwrite=True)