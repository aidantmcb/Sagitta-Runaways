import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt

#### Load files for clustering and age predictions
kc19_fname = 'c:/users/sahal/desktop/ML Astro/Data/final6age.fits'
sagitta_fname = 'c:/users/sahal/desktop/sagitta-edr3.fits'

def maketable(fname):
    data = np.array(fits.open(fname)[1].data)
    return pd.DataFrame(data.byteswap().newbyteorder())

#### Find average parameters for individual star forming regions

#Takes kc19 as table
def clusterparams(table):
    avg_plx = np.mean(table['PARALLAX'])
    avg_l, avg_b = (np.mean(table['L']), np.mean(table['B']))
    std_l, std_b = (np.std(table['L'], ddof=1), np.std(table['B'], ddof =1))
    return avg_plx, avg_l, avg_b, std_l, std_b

def clustermotionparams(table):
    avg_pml, avg_pmb = (np.mean(table['VLSRL']), np.mean(table['VLSRB']))
    std_pml, std_pmb = (np.std(table['VLSRL'], ddof = 1), np.std(table['VLSRB'], ddof = 1))
    return avg_pml, avg_pmb, std_pml, std_pmb

def getregion(table, spatial, motion):
    avg_plx, avg_l, avg_b, std_l, std_b = spatial
    avg_pml, avg_pmb, std_pml, std_pmb = motion
    rad = np.sqrt(std_l**2 + std_b**2)
    circ = np.where((np.sqrt((avg_l-table['l'])**2 + (avg_b-table['b'])**2)) < rad)[0]

    return circ

#### Establish a probability of individual stars within some region of being cluster members
def __main__():
    kc19 = maketable(kc19_fname)
    sagitta = maketable(sagitta_fname)
    ids = [13]
    for id in ids:
        cluster = kc19.iloc[np.where(kc19['ID']== id)[0]]
        spatialparams = clusterparams(cluster)
        motionparams = clustermotionparams(cluster)
        region = getregion(sagitta, spatialparams, motionparams)
        k = sagitta.iloc[region]
        plt.scatter(k['l'], k['b'], s=0.1)
        plt.show()

__main__()