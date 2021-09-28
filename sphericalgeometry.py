import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from astropy.table import Table

from tol_colors import tol_cmap
from RunawayPlots import plotSky
from RunawayPlots import plotSTILTS

kc19_fname = "c:/users/sahal/desktop/ML Astro/Data/final6age.fits"
sagitta_fname = "c:/users/sahal/desktop/Sagitta_HDBSCAN_named_.fits"

def maketable(fname):
    data = np.array(fits.open(fname)[1].data)
    table = pd.DataFrame(data.byteswap().newbyteorder())
    for col in table:
        table.rename(columns={col: col.lower()}, inplace=True)
    if 'pms' in table.columns:
        table = table.rename(columns = {'pms': 'yso'})
    if 'clusterer' in table.columns:
        table = table.rename(columns = {'clusterer': 'labels'})
    # if 'pms_mean' in table.columns:
    #     table.rename(columns = {'pms_mean': 'yso'})
    return table

def motions(pml, pmb, l, b):
    

