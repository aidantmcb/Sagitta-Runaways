import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.table import Table
import hdbscan

fname = 'C:/users/sahal/Desktop/Sagitta-edr3_HDBSCAN_.fits'

def maketable(fname):
    data = np.array(fits.open(fname)[1].data)
    table = pd.DataFrame(data.byteswap().newbyteorder())
    for col in table:
        table.rename(columns={col: col.lower()}, inplace=True)
    return table

df = maketable(fname)
df = df.iloc[np.where(df['parallax']>1)[0]]
velo_l = df['vlsrl'] / 13600000 * np.pi/180 * 1000 / df['parallax'] * (3.086e13) / (3.154e7) #km/s
velo_b = df['vlsrb'] / 13600000 * np.pi/180 * 1000 / df['parallax'] * (3.086e13) / (3.154e7) #km/s
df['velocity_l'], df['velocity_b'] = velo_l, velo_b
df['dist'] = 1000/df['parallax']

clusterer = hdbscan.HDBSCAN()
clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']])
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_1'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN()
clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb', 'age']])
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_2'] = clusterer.labels_

# clusterer = hdbscan.HDBSCAN(min_cluster_size = 10)
# clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']])
# print('Labels:', np.max(clusterer.labels_) + 1)
# df['labels_3'] = clusterer.labels_

# clusterer = hdbscan.HDBSCAN(min_cluster_size = 10)
# clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb', 'age']])
# print('Labels:', np.max(clusterer.labels_) + 1)
# df['labels_4'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN()
clusterer.fit(df[['l', 'b', 'parallax', 'velocity_l', 'velocity_b']])
print('Labels:', np.max(clusterer.labels_) + 1)
df['labels_5'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN(cluster_selection_method = 'leaf')
clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']])
print('Labels:', np.max(clusterer.labels_) + 1)
df['labels_6'] = clusterer.labels_

# clusterer = hdbscan.HDBSCAN(cluster_selection_method = 'leaf')
# clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb', 'age']])
# print('Labels:', np.max(clusterer.labels_) + 1)
# df['labels_7'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN()
lowest = np.min(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']]) 
norm = np.max(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']]) - np.min(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']]) 
data = (df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']]-lowest)/norm
clusterer.fit(data)
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_8'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN(min_cluster_size = 20, min_samples = 2)
clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']])
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_9'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN()
clusterer.fit(df[['l', 'b', 'dist', 'vlsrl', 'vlsrb']])
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_10'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN(alpha = 0.5)
clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']])
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_11'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN(min_cluster_size=3)
clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']])
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_12'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN(min_cluster_size = 3)
lowest = np.min(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']]) 
norm = np.max(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']]) - np.min(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']]) 
data = (df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']]-lowest)/norm
clusterer.fit(data)
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_13'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN(min_cluster_size=3)
clusterer.fit(df[['l', 'b', 'vlsrl', 'vlsrb']])
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_14'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN(min_cluster_size=3, alpha=0.6)
clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']])
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_15'] = clusterer.labels_

clusterer = hdbscan.HDBSCAN(min_cluster_size=2)
clusterer.fit(df[['l', 'b', 'parallax', 'vlsrl', 'vlsrb']])
print('Labels:', np.max(clusterer.labels_)+ 1)
df['labels_16'] = clusterer.labels_

for i in range(1, 17):
    labelsname = 'labels_' + str(i)
    colname = 'c' + str(i)
    df[colname] = np.where(df[labelsname]!=-1, True, False)

df = df.iloc[np.where(df['yso']>0.85)[0]]

tab = Table.from_pandas(df)
tab.write('c:/users/sahal/Desktop/Sagitta_HDBSCAN_.fits',overwrite=True)


