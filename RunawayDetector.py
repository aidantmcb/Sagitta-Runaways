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
sagitta_fname = "c:/users/sahal/desktop/Sagitta_HDBSCAN_named.fits"


def maketable(fname):
    data = np.array(fits.open(fname)[1].data)
    table = pd.DataFrame(data.byteswap().newbyteorder())
    for col in table:
        table.rename(columns={col: col.lower()}, inplace=True)
    if 'pms' in table.columns:
        table.rename(columns = {'pms': 'yso'})
    return table


#### Find average parameters for individual star forming regions

# Takes kc19 as table
def clusterspatialparams(table):
    avg_plx = np.mean(table["parallax"])
    avg_l, avg_b = (np.mean(table["l"]), np.mean(table["b"]))
    std_l, std_b = (np.std(table["l"], ddof=1), np.std(table["b"], ddof=1))
    avg_age = np.mean(table['age'])
    return avg_plx, avg_l, avg_b, std_l, std_b, avg_age


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

def get_dist(l, b, avg_l, avg_b):
        d = np.sqrt((l-avg_l)**2 + (b-avg_b)**2)
        return d

def get_relative_velocity(pm_b, pm_l, parallax, avg_pml, avg_pmb):
    pm_l_c, pm_b_c = pm_l - avg_pml, pm_b - avg_pmb
    pm_mag = np.sqrt(pm_l_c**2 + pm_b_c**2) / 13600000 * np.pi/180
    v_relative = pm_mag * 1000 / parallax * (3.086e13) / (3.154e7) #km/s
    return v_relative

def get_tracebacktime(l, b, pm_l, pm_b, avg_l, avg_b, avg_pml, avg_pmb):
    pm_mag = np.sqrt((pm_l - avg_pml)**2 + (pm_b - avg_pmb)**2)
    tbt = np.sqrt((l-avg_l)**2 + (b-avg_b)**2) / pm_mag * 13600000
    return tbt


#### Establish a probability of individual stars within some region of being cluster members
def main():
    # kc19 = maketable(kc19_fname) #Get table for cluster averages
    # kc19 = kc19.iloc[np.where(kc19['parallax'] > 1)[0]]
    # young_clusters = np.where(np.power(10,(kc19['age']-6))<45)[0]
    # labels = np.sort(np.unique(kc19['labels'].iloc[young_clusters]))
    # print(labels)
    # print('Regions:', len(labels))

    sagitta = maketable(sagitta_fname) #Get table for model outputs    
    sagitta = sagitta.iloc[
        np.where((sagitta["yso"] > 0.90) & (sagitta['parallax'] > 1))[0]
    ]
    if 'l1' in sagitta.columns:
        sagitta.drop(columns=['l1']) #reset this column if it already exists
    sagitta['l1'] = sagitta['l']
    over180 = np.where(sagitta['l1']>180)[0]
    sagitta['l1'].iloc[over180] = sagitta['l1'].iloc[over180]-360


    labels = np.sort(np.unique(sagitta['labels']))
    if -1 in labels:
        labels = np.delete(labels, 0)
    print(len(labels))
    
    addcols = ['traceback','avg_l', 'avg_b', 'avg_pml', 'avg_pmb', 'avg_age', 'cluster_id', 'cluster_label']
    output_frame = pd.DataFrame(columns=(list(sagitta.columns)+addcols))
    # Orion: 13, Serpens: 23 - labels: Orion 29, Oph 126,

    for label in labels:
        print(label)
        # cluster = kc19.iloc[np.where(kc19["labels"] == label)[0]]
        cluster = sagitta.iloc[np.where(sagitta["labels"] == label)[0]]

        avg_plx, avg_l, avg_b, std_l, std_b, avg_age = clusterspatialparams(cluster)
        avg_pml, avg_pmb, std_pml, std_pmb = clustermotionparams(cluster)

        if 'id' in cluster.columns: 
            ids = np.min(np.unique(cluster['id']))
        else: 
            ids = 00

        l = 'l'
        avg_l_save = avg_l

        if avg_l < 90 or 360 - avg_l < 90:
            l = 'l1'
            if 360 - avg_l < 90:
                avg_l = avg_l - 360 #Needed
 
        # print(np.power(10,sagitta['age'].iloc[5]-6), np.power(10, avg_age))

        # break
        fn = lambda a : -(1/(7.7-6))*(a-6)+2 #Double age range for younger ages, but shrink that for older ones
        tab = sagitta.iloc[np.where(np.power(10, sagitta['age']-6) < np.power(10, avg_age-6)*fn(avg_age))[0]]

        #NEW:
        tab = tab.iloc[np.where((tab['labels'] == label) | (tab['labels'] == -1))[0]]

        alignment = getalignment(
            tab[l],
            tab["b"],
            tab["vlsrl"],
            tab["vlsrb"],
            avg_l,
            avg_b,
            avg_pml,
            avg_pmb,
        )

        max_offset = 0.1
        align_fn = lambda d : 1.025**-d
        aligned = np.where(alignment > 1-(max_offset * align_fn(get_dist(tab[l],tab['b'],avg_l, avg_b))))[0]

        

        tab = tab.iloc[aligned] #SELECTION 1: alignment > 1-max_offset

        
        velos = get_relative_velocity(tab['vlsrl'], tab['vlsrb'], tab['parallax'], avg_pml, avg_pmb)
        tab = tab.iloc[np.where(velos>10)[0]] # SELECTION  2: velocities > 10 km/s

        tbt = get_tracebacktime(tab[l], tab['b'], tab['vlsrl'],  tab['vlsrb'], avg_l, avg_b, avg_pml, avg_pmb) 
        tab['traceback'] = np.log10(tbt)
        tab = tab.iloc[np.where(tbt <  np.power(10, avg_age))[0]]

        #tab = tab.iloc[np.where(np.abs(tab['parallax'] - avg_plx) < tab['eparallax'])]

        # plotSTILTS(avg_l , avg_b, avg_pml, avg_pmb, .80, 6.8, align_threshold = max_offset, aux_max = 1e7, l = l)

        tab['avg_l'] = np.float(avg_l_save)
        tab['avg_b'] = np.float(avg_b)
        tab['avg_pml'] = np.float(avg_pml)
        tab['avg_pmb'] = np.float(avg_pmb)
        tab['avg_age'] = np.float(avg_age)
        tab['cluster_id'] = int(ids)
        tab['cluster_label'] = int(label)
        output_frame = output_frame.append(tab)

    for column in output_frame.columns:
        if (output_frame[column].dtype == 'object') & (column not in ['name', 'plotname']):
            output_frame[column] = output_frame[column].astype('float')
    
    t = Table.from_pandas(output_frame)
    t.write('c:/users/sahal/desktop/RunawayDetector_4-15a.fits', overwrite=True)


if __name__ == "__main__":
    main()
