import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.table import Table
import hdbscan

fname = 'C:/users/sahal/Desktop/Sagitta-edr3-final6.fits'

def maketable(fname):
    data = np.array(fits.open(fname)[1].data)
    table = pd.DataFrame(data.byteswap().newbyteorder())
    for col in table:
        table.rename(columns={col: col.lower()}, inplace=True)
    if 'pms' in table.columns:
        table.rename(columns = {'pms': 'yso'}, inplace=True)
    return table

df = maketable(fname)
df = df.iloc[np.where(df['parallax']>1)[0]]
# velo_l = df['vlsrl'] / 13600000 * np.pi/180 * 1000 / df['parallax'] * (3.086e13) / (3.154e7) #km/s
# velo_b = df['vlsrb'] / 13600000 * np.pi/180 * 1000 / df['parallax'] * (3.086e13) / (3.154e7) #km/s
# df['velocity_l'], df['velocity_b'] = velo_l, velo_b
# df['dist'] = 1000/df['parallax']



clusterer1 = hdbscan.HDBSCAN(min_cluster_size = 20)
clusterer1.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']])
print('Labels:', np.max(clusterer1.labels_) + 1)
df['labels'] = clusterer1.labels_

clusterer2 = hdbscan.HDBSCAN(min_cluster_size = 10)
unclustered = np.where((df['labels'] == -1) & (df['yso']>.9))[0]
clusterer2.fit(df[['l','b','parallax','vlsrl','vlsrb']].iloc[unclustered])
print('Labels:', np.max(clusterer2.labels_) + 1)
labels2 = clusterer2.labels_
labels2[np.where(labels2 != -1)[0]] = labels2[np.where(labels2 != -1)[0]] + np.max(df['labels']) + 1
df['labels'].iloc[unclustered] = labels2

tab = Table.from_pandas(df)
tab.write('c:/users/sahal/Desktop/Sagitta_HDBSCAN-final6.fits', overwrite=True)
