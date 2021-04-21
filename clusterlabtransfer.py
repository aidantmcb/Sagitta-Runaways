import numpy as np
import pandas as pd
from astropy.io import fits
import os
from astropy.table import Table

fname = 'c:/Users/sahal/Desktop/Sagitta_HDBSCAN-final6.fits'
data = np.array(fits.open(fname)[1].data)
df = pd.DataFrame(data.byteswap().newbyteorder())

# f6fname = 'c:/users/sahal/Desktop/final6age.fits'
# f6data = np.array(fits.open(fname)[1].data)
# f6 = pd.DataFrame(data.byteswap().newbyteorder())

df['f6_name'] = np.zeros(len(df))

for label in np.unique(df['labels']):
    print(label)
    x = np.where(df['labels']==label)[0]
    cluster = df.iloc[x]
    counts = cluster['name'].iloc[np.where(cluster['name'] != b' ')[0]].value_counts()

    name = b' '
    if counts.size > 0:
        name = counts.index[0]
    df['f6_name'].iloc[x] = name

# already_named = np.where(df['name'] != b' ')
# df['f6_name'].iloc[already_named] = df['name'].iloc[already_named]

tab = Table.from_pandas(df)
tab.write('c:/users/sahal/desktop/Sagitta_HDBSCAN_named.fits', overwrite = True)