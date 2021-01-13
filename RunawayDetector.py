import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from tol_colors import tol_cmap
from astropy.table import Table

rainbow = tol_cmap("rainbow_PuRd")
rainbow_r = rainbow.reversed()

#### Load files for clustering and age predictions
kc19_fname = "c:/users/sahal/desktop/ML Astro/Data/final6age.fits"
sagitta_fname = "c:/users/sahal/desktop/sagitta-edr3.fits"


def maketable(fname):
    data = np.array(fits.open(fname)[1].data)
    table = pd.DataFrame(data.byteswap().newbyteorder())
    for col in table:
        table.rename(columns={col: col.lower()}, inplace=True)
    return table

#### Find average parameters for individual star forming regions

# Takes kc19 as table
def clusterparams(table):
    avg_plx = np.mean(table["parallax"])
    avg_l, avg_b = (np.mean(table["l"]), np.mean(table["b"]))
    std_l, std_b = (np.std(table["l"], ddof=1), np.std(table["b"], ddof=1))
    return avg_plx, avg_l, avg_b, std_l, std_b


def clustermotionparams(table):
    avg_pml, avg_pmb = (np.mean(table["vlsrl"]), np.mean(table["vlsrb"]))
    std_pml, std_pmb = (np.std(table["vlsrl"], ddof=1), np.std(table["vlsrb"], ddof=1))
    return avg_pml, avg_pmb, std_pml, std_pmb


def getregion(table, spatial, motion, radius_modifier=1):
    avg_plx, avg_l, avg_b, std_l, std_b = spatial
    avg_pml, avg_pmb, std_pml, std_pmb = motion
    rad = np.sqrt(std_l ** 2 + std_b ** 2)
    circ = np.where(
        (np.sqrt((avg_l - table["l"]) ** 2 + (avg_b - table["b"]) ** 2)) < rad * radius_modifier
    )[0]
    return circ


#### Establish a probability of individual stars within some region of being cluster members
def main():
    kc19 = maketable(kc19_fname)
    sagitta = maketable(sagitta_fname)
    sagitta = sagitta.iloc[np.where((sagitta['yso']>.80) & (np.power(10,sagitta['age']-6)<45))[0]]
    ids = [13] #Orion: 13, Serpens: 23 - 
    for id in ids:
        cluster = kc19.iloc[np.where(kc19["labels"] == 29)[0]]
        spatialparams = clusterparams(cluster)
        motionparams = clustermotionparams(cluster)
        region_indices = getregion(sagitta, spatialparams, motionparams, radius_modifier = 4)
        region = sagitta.iloc[region_indices]

        fig = plt.figure()
        ax = plt.subplot()
        points = ax.scatter(region["l"], region["b"], c=np.power(10, region['age']-6), s=0.5, cmap=rainbow_r)
        fig.colorbar(points, ax = ax, label = "Age (Myr)")
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(xmax, xmin)
        plt.show()

        fig = plt.figure()
        ax = plt.subplot()
        avg_pml, avg_pmb, std_pml, std_pmb = motionparams[0], motionparams[1], motionparams[2], motionparams[3]
        ax.scatter(region['vlsrl']-avg_pml, region['vlsrb']-avg_pmb, s=0.8)
        k=3
        ax.plot([k*std_pml, k*std_pml, k*-std_pml, k*-std_pml, k*std_pml], 
        [k*std_pmb, k*-std_pmb, k*-std_pmb, k*std_pmb, k*std_pmb],color='orange')
        ax.set_xlabel("VLSRL")
        ax.set_ylabel("VLSRB")

        x = np.where((np.abs(region['vlsrl']-avg_pml)>k*std_pml) | (np.abs(region['vlsrb']-avg_pmb)> k*std_pmb))[0] 
        print(len(x))
        OUT = region.iloc[x]
        ax.scatter(OUT['vlsrl']-avg_pml, OUT['vlsrb']-avg_pmb, s=0.8, c='r')
        plt.show()
        

        t = Table.from_pandas(OUT)
        t.write('c:/users/sahal/desktop/OUT.fits', overwrite=True)




if __name__ == "__main__":
    main()