import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import LSR, ICRS, SkyCoord
from astropy import units as u

oriontab = '/Users/aidanmcbride/Documents/Sagitta-Runaways/master_revised.fits'
runawaytab = '/Users/aidanmcbride/Documents/Sagitta-Runaways/SphericalMotions/muCol_aeAur.fits'

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

ori = maketable(oriontab)
run = maketable(runawaytab)

# c1_l = SkyCoord(ra=run['ra'].iloc[0]*u.degree, dec=run['dec'].iloc[0]*u.degree, frame='icrs',pm_ra_cosdec=run['pmra'].iloc[0]*u.mas/u.yr, pm_dec=run['pmdec'].iloc[0]*u.mas/u.yr,distance=1000/run['parallax'].iloc[0]*u.pc,radial_velocity=0*u.km/u.s)
# c2_l = SkyCoord(ra=run['ra'].iloc[1]*u.degree, dec=run['dec'].iloc[1]*u.degree, frame='icrs',pm_ra_cosdec=run['pmra'].iloc[1]*u.mas/u.yr, pm_dec=run['pmdec'].iloc[1]*u.mas/u.yr,distance=1000/run['parallax'].iloc[1]*u.pc,radial_velocity=0*u.km/u.s)
# c1 = c1_l.transform_to('lsr')
# c2 = c2_l.transform_to('lsr')

ae_c = SkyCoord(ra=79.0756206*u.degree, dec=34.31231862*u.degree, frame='icrs',pm_ra_cosdec=-4.440*u.mas/u.yr, pm_dec=43.368*u.mas/u.yr,distance=1000/2.4642*u.pc,radial_velocity=0*u.km/u.s)
ae_l=ae_c.transform_to(LSR())
mu_c = SkyCoord(ra=86.49956401*u.degree, dec=-32.30643465*u.degree, frame='icrs',pm_ra_cosdec=2.988*u.mas/u.yr, pm_dec=-22.030*u.mas/u.yr,distance=1000/2.1476*u.pc,radial_velocity=0*u.km/u.s)
mu_l=mu_c.transform_to(LSR())

c1 = ae_l
c2 = mu_l

cluster = SkyCoord(ra=83.8*u.degree, dec=-5.4*u.degree, frame='icrs',pm_ra_cosdec=0.8*u.mas/u.yr, pm_dec=2.8*u.mas/u.yr,distance=390*u.pc,radial_velocity=0*u.km/u.s)

def bearing(lat1, lon1, lat2, lon2):
    lat1 = lat1 * np.pi / 180
    lon1 = lon1 * np.pi / 180
    lat2 = lat2 * np.pi / 180
    lon2 = lon2 * np.pi / 180
    y = np.sin(lon2-lon1) * np.cos(lat2)
    x = np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)
    # return np.arctan2(y,x)
    return np.arctan2(x,y)

lon1 = c1.ra.value #* np.pi/180
lat1 = c1.dec.value #* np.pi/180
lon2 = cluster.ra.value #* np.pi/180
lat2 = cluster.dec.value #* np.pi/180
lon3 = (c1.ra.value - (c1.pm_ra_cosdec.value - cluster.pm_ra_cosdec.value)/np.cos(c1.dec.value * np.pi/180)/1e3/3600) #* np.pi/180
lat3 = (c1.dec.value - (c1.pm_dec.value - cluster.pm_dec.value)/1e3/3600) #* np.pi / 180

b1 = bearing(lat1, lon1, lat2, lon2)
b2 = bearing(lat1, lon1, lat3, lon3)
print(b1)
print(b2)
print((b1-b2)*180/np.pi)
print(np.cos((b1-b2))) #* 180 / np.pi))