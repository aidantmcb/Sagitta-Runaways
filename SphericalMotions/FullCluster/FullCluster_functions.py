import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

from astropy.table import Table
from astropy.io import fits
from astropy.coordinates import LSR, ICRS, SkyCoord
from astropy import units as u

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

def within_distance(tab, avg_plx, thrshld = 100): #requires within 100 pc
    dist = 1000/tab['parallax']
    avg_dist = 1000/avg_plx
    keep = np.where(np.abs(dist-avg_dist)<thrshld)[0]
    return tab.iloc[keep]

#####
#Uses spherical geometry calculations with math found at http://www.movable-type.co.uk/scripts/latlong.html
#####

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