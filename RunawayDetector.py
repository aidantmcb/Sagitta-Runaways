import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from astropy.table import Table

from tol_colors import tol_cmap
from RunawayPlots import plotSky
from RunawayPlots import plotSTILTS

rainbow = tol_cmap("rainbow_PuRd")
rainbow_r = rainbow.reversed()

#### Load files for clustering and age predictions
kc19_fname = "c:/users/sahal/desktop/ML Astro/Data/final6age.fits"
sagitta_fname = "c:/users/sahal/desktop/sagitta_edr3-l1.fits"


def maketable(fname):
    data = np.array(fits.open(fname)[1].data)
    table = pd.DataFrame(data.byteswap().newbyteorder())
    for col in table:
        table.rename(columns={col: col.lower()}, inplace=True)
    return table


#### Find average parameters for individual star forming regions

# Takes kc19 as table
def clusterspatialparams(table):
    avg_plx = np.mean(table["parallax"])
    avg_l, avg_b = (np.mean(table["l"]), np.mean(table["b"]))
    std_l, std_b = (np.std(table["l"], ddof=1), np.std(table["b"], ddof=1))
    return avg_plx, avg_l, avg_b, std_l, std_b


def clustermotionparams(table):
    avg_pml, avg_pmb = (np.mean(table["vlsrl"]), np.mean(table["vlsrb"]))
    std_pml, std_pmb = (np.std(table["vlsrl"], ddof=1), np.std(table["vlsrb"], ddof=1))
    return avg_pml, avg_pmb, std_pml, std_pmb


def getalignment(l, b, pml, pmb, avg_l, avg_b, avg_pml, avg_pmb):
    angle = np.arctan2(pml - avg_pml, pmb - avg_pmb) - np.arctan2(l - avg_l, b - avg_b)
    return np.cos(angle)

# atan2(vlsrl + 3.4, vlsrb + 4.6) - atan2(l1 +6.8, b - 17.25)
def getregion(table, spatial, motion, radius_modifier=1):
    avg_plx, avg_l, avg_b, std_l, std_b = spatial
    avg_pml, avg_pmb, std_pml, std_pmb = motion
    rad = np.sqrt(std_l ** 2 + std_b ** 2)
    circ = np.where(
        (np.sqrt((avg_l - table["l"]) ** 2 + (avg_b - table["b"]) ** 2))
        < rad * radius_modifier
    )[0]
    return circ

def get_relative_velocity(pm_b, pm_l, parallax, avg_pml, avg_pmb):
    pm_l_c, pm_b_c = pm_l - avg_pml, pm_b - avg_pmb
    pm_mag = np.sqrt(pm_l_c**2 + pm_b_c**2) / 13600000 * np.pi/180
    v_relative = pm_mag * 1000 / parallax * (3.086e13) / (3.154e7) #km/s
    return v_relative

#### Establish a probability of individual stars within some region of being cluster members
def main():
    kc19 = maketable(kc19_fname)
    sagitta = maketable(sagitta_fname)
    sagitta = sagitta.iloc[
        np.where((sagitta["yso"] > 0.80) & (np.power(10, sagitta["age"] - 6) < 45))[0]
    ]
    ids = [13]  # Orion: 13, Serpens: 23 - labels: Orion 29, Oph 126
    for id in ids:
        cluster = kc19.iloc[np.where(kc19["labels"] == 126)[0]]
        avg_plx, avg_l, avg_b, std_l, std_b = clusterspatialparams(cluster)
        avg_pml, avg_pmb, std_pml, std_pmb = clustermotionparams(cluster)

        l = 'l'
        if avg_l < 90 or 360 - avg_l < 90:
            l = 'l1'
            if 360 - avg_l < 90:
                avg_l = avg_l - 360 #Needed


        alignment = getalignment(
            sagitta[l],
            sagitta["b"],
            sagitta["vlsrl"],
            sagitta["vlsrb"],
            avg_l,
            avg_b,
            avg_pml,
            avg_pmb,
        )
        max_offset = 0.1
        aligned = np.where(alignment > 1-max_offset)[0]
        # print('Aligned Stars:', len(aligned))

        plotSTILTS(avg_l , avg_b, avg_pml, avg_pmb, .80, 6.8, radius = 35, align_threshold = max_offset, l = l)


        # tab = sagitta.iloc[aligned]
        # plotSky(tab, alignment[aligned], avg_l, avg_b)
        
        # velos = get_relative_velocity(tab['vlsrl'], tab['vlsrb'], tab['parallax'], avg_pml, avg_pmb)
        # plotSky(tab, np.log10(velos), avg_l, avg_b)

        # t = Table.from_pandas(tab)
        # t.write('c:/users/sahal/desktop/OUT.fits', overwrite=True)


if __name__ == "__main__":
    main()
