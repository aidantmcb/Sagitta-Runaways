import numpy as np
import pandas as pd
from astropy.io import fits
import os
from astropy.table import Table

fname = 'c:/Users/sahal/Desktop/RunawayDetector_4-27TEST.fits'
data = np.array(fits.open(fname)[1].data)
df = pd.DataFrame(data.byteswap().newbyteorder())
print('Unique Sources: ' + str(len(np.unique(df['source_id'], return_index = True)[0])))

def get_dist(l, b, avg_l, avg_b):
        d = np.sqrt((l-avg_l)**2 + (b-avg_b)**2)
        return d

# def getl1(tabl):
#         over180 = np.where(tabl>180)[0]
#         tabl.iloc[over180] = tabl.iloc[over180]-360
#         return tabl

def getl1(tabl):
        avgl = np.array(tabl)
        i = np.where(360 - avgl < 90)[0]
        avgl[i] = avgl[i] - 360
        return avgl 
        # if 360 - avg_l < 90:
        #         avg_l = avg_l - 360 #Needed

l_projected = df['l'] - (df['vlsrl'] - df['avg_pml'] / 13600000 * np.power(10,df['traceback']))
cond = np.where((df['avg_l'] < 90) | (360 - df['avg_l'] < 90 ))[0]
l_projected.iloc[cond] = df['l1'].iloc[cond] - (df['vlsrl'].iloc[cond]-df['avg_pml'].iloc[cond]) / 13600000 * np.power(10, df['traceback'].iloc[cond])
b_projected = df['b'] - (df['vlsrb']-df['avg_pmb']) / 13600000 * np.power(10, df['traceback'])

print(len(cond))
#Select by stdev of l and b - only possible for tables made after 4-20-21
sigma = 2
within_l = np.where(np.abs(l_projected - df['avg_l']) < sigma * df['std_l'])[0]
within_l1 = np.where(np.abs(l_projected.iloc[cond] - getl1(df['avg_l'].iloc[cond])) < sigma * df['std_l'].iloc[cond])[0]
within_l_l1 = np.union1d(within_l, within_l1)
within_b = np.where(np.abs(b_projected - df['avg_b']) < sigma * df['std_b'])[0]
retain = np.intersect1d(within_l_l1, within_b)
df = df.iloc[retain]


#Select by closest distance to center of cluster
# dist = np.array(get_dist(l_projected, b_projected, getl1(df['avg_l']), df['avg_b']))
# df = df.iloc[np.argsort(dist)]


#Select by shortest traceback time
df = df.iloc[np.argsort(df['traceback'])]

#Select by closest traceback to age
# timediff = np.abs(df['age']-df['traceback'])
# df = df.iloc[np.argsort(timediff)]


x = np.unique(df['source_id'], return_index=True)[1]
df = df.iloc[x]
df = df.iloc[np.argsort(df['cluster_label'])]
print('Final length: ' + str(len(df)))

t = Table.from_pandas(df) 
t.write('c:/users/sahal/desktop/RunawayDetector_4-27TEST-tracebacksorted.fits', overwrite=True)