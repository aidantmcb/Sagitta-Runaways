import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import LSR, ICRS, SkyCoord
from astropy import units as u

from RunawayDetectorspherical import (maketable, clusterspatialparams, clustermotionparams, getalignment, getregion, get_dist,
                get_relative_velocity, get_tracebacktime, within_distance)

#####
#Uses spherical geometry calculations with math found at http://www.movable-type.co.uk/scripts/latlong.html
#####

kc19_fname = '/Users/aidanmcbride/Documents/Sagitta-Runaways/final6age.fits'
sagitta_fname = '/Users/aidanmcbride/Documents/Sagitta-Runaways/Sagitta_HDBSCAN_named_.fits'


def bearing(lat1, lon1, lat2, lon2):
    lat1 = lat1 * np.pi / 180
    lon1 = lon1 * np.pi / 180
    lat2 = lat2 * np.pi / 180
    lon2 = lon2 * np.pi / 180
    y = np.sin(lon2-lon1) * np.cos(lat2)
    x = np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)
    # return np.arctan2(y,x)
    return np.arctan2(x,y)

def spheredist(lat1, lon1, lat2, lon2):
    lat1 = lat1 * np.pi / 180
    lon1 = lon1 * np.pi / 180
    lat2 = lat2 * np.pi / 180
    lon2 = lon2 * np.pi / 180
    dlat = lat2 - lat1 
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1)*np.cos(lat2)*np.sin(dlon/2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
    return c * 180 / np.pi

def sphere_tracebacktime(dist, pm_l, pm_b, avg_pml, avg_pmb):
    pm_mag = np.sqrt((pm_l - avg_pml)**2 + (pm_b - avg_pmb)**2)
    tbt = dist / pm_mag * 13600000
    return tbt


def main():
    sagitta = maketable(sagitta_fname) #Get table for model outputs    
    sagitta.rename(columns = {'pms': 'yso'})
    sagitta = sagitta.iloc[
        np.where((sagitta['yso'] > 0.90) & (sagitta['parallax'] > 1))[0]
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

    addcols = ['traceback','avg_l', 'avg_b', 'std_l', 'std_b' 'avg_pml', 'avg_pmb','avg_plx', 
                    'avg_age', 'cluster_id', 'cluster_label', 'alignment', 'gcdistance', 'olddist']
    output_frame = pd.DataFrame(columns=(list(sagitta.columns)+addcols))

    for label in labels: 
        print(label)
        cluster = sagitta.iloc[np.where(sagitta["labels"] == label)[0]]

        avg_plx, avg_l, avg_b, std_l, std_b, avg_age = clusterspatialparams(cluster)
        avg_pml, avg_pmb, std_pml, std_pmb = clustermotionparams(cluster)

        if 'id' in cluster.columns: 
            ids = np.min(np.unique(cluster['id']))
        else: 
            ids = 00

        # sagitta['l'] = sagitta['l'] + 20 #added to test chance alignments

        l = 'l'
        avg_l_save = avg_l

        if avg_l < 90 or 360 - avg_l < 90:
            l = 'l1'
            if 360 - avg_l < 90:
                avg_l = avg_l - 360 #Needed

        fn = lambda a : -(1/(7.7-6))*(a-6)+2 #Double age range for younger ages, but shrink that for older ones
        tab = sagitta.iloc[np.where(np.power(10, sagitta['age']-6) < np.power(10, avg_age-6)*fn(avg_age))[0]]
        tab = tab.iloc[np.where((tab['labels'] == label) | (tab['labels'] == -1))[0]]

        lon3=(tab[l]-(tab['vlsrl']-avg_pml)/1e3/3600)#*np.pi/180
        # lon3=(tab['l']-(tab['vlsrl']-avg_pml)/1e3/3600)#*np.pi/180

        lat3=(tab['b']-(tab['vlsrb']-avg_pmb)/1e3/3600)#*np.pi/180

        # bearing_to = bearing(tab[l], tab['b'], avg_l, avg_b)
        # bearing_of = bearing(lat3, lon3, avg_l, avg_b)

        bearing_to = bearing(tab['b'], tab[l], avg_b, avg_l)
        bearing_of = bearing(tab['b'], tab[l], lat3, lon3)
        alignment = np.cos(bearing_of - bearing_to)

        max_offset = 0.1
        align_fn = lambda d : 1.025**-d
        aligned = np.where(alignment > 1-(max_offset * align_fn(get_dist(tab[l],tab['b'],avg_l, avg_b))))[0] 
        tab['alignment'] = alignment

        # tab['gcdist'] = spheredist(tab['b'], tab[l], avg_b, avg_l)#added 11-4-21
        # tab['olddist'] = get_dist(tab[l],tab['b'], avg_l, avg_b)
        # tab = tab.iloc[aligned] #SELECTION 1: alignment > 1-max_offset
    
        velos = get_relative_velocity(tab['vlsrl'], tab['vlsrb'], tab['parallax'], avg_pml, avg_pmb)
        tab = tab.iloc[np.where(velos>10)[0]] # SELECTION  2: velocities > 10 km/s

        tbt = get_tracebacktime(tab[l], tab['b'], tab['vlsrl'],  tab['vlsrb'], avg_l, avg_b, avg_pml, avg_pmb) 
        # tbt = sphere_tracebacktime(tab['gcdist'], tab['vlsrl'], tab['vlsrb'], avg_pml, avg_pmb)
        tab['traceback'] = np.log10(tbt)
        tab = tab.iloc[np.where(tbt <  np.power(10, avg_age))[0]]
        tab = within_distance(tab, avg_plx, thrshld=100)

        #tab = tab.iloc[np.where(np.abs(tab['parallax'] - avg_plx) < tab['eparallax'])]

        # plotSTILTS(avg_l , avg_b, avg_pml, avg_pmb, .80, 6.8, align_threshold = max_offset, aux_max = 1e7, l = l)

        tab['avg_l'] = np.float(avg_l_save)
        tab['avg_b'] = np.float(avg_b)
        tab['std_l'] = np.float(std_l)
        tab['std_b'] = np.float(std_b)
        tab['avg_pml'] = np.float(avg_pml)
        tab['avg_pmb'] = np.float(avg_pmb)
        tab['avg_plx'] = np.float(avg_plx)
        tab['avg_age'] = np.float(avg_age)
        tab['cluster_id'] = int(ids)
        tab['cluster_label'] = int(label)
        output_frame = output_frame.append(tab)

    for column in output_frame.columns:
        if (output_frame[column].dtype == 'object') & (column not in ['name', 'f6_name' ,'plotname']):
            output_frame[column] = output_frame[column].astype('float')
    
    t = Table.from_pandas(output_frame)
    t.write('/Users/aidanmcbride/Documents/Sagitta-Runaways/11-5-21_CHANCE.fits', overwrite=True)
if __name__ == '__main__':
    main()