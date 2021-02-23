import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.table import Table
import hdbscan

fname = 'C:/users/sahal/Desktop/Sagitta-edr3.fits'

def maketable(fname):
    data = np.array(fits.open(fname)[1].data)
    table = pd.DataFrame(data.byteswap().newbyteorder())
    for col in table:
        table.rename(columns={col: col.lower()}, inplace=True)
    return table

df = maketable(fname)
df = df.iloc[np.where(df['parallax']>1)[0]]
# velo_l = df['vlsrl'] / 13600000 * np.pi/180 * 1000 / df['parallax'] * (3.086e13) / (3.154e7) #km/s
# velo_b = df['vlsrb'] / 13600000 * np.pi/180 * 1000 / df['parallax'] * (3.086e13) / (3.154e7) #km/s
# df['velocity_l'], df['velocity_b'] = velo_l, velo_b
# df['dist'] = 1000/df['parallax']

