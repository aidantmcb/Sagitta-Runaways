import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import LSR, ICRS, SkyCoord
from astropy import units as u


oriontab = 'c:/users/sahal/documents/OrionProperMotions/master_revised.fits'
runawaytab = 'C:/Users/sahal/Documents/ClusterKinematics/SphericalMotions/muCol_aeAur.fits'

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
# print(run.columns)

avgra, avgdec = (np.mean(ori['ra']), np.mean(ori['dec']))


d13 = np.sqrt((run['ra']-avgra)**2 + (run['dec']-avgdec)**2)  * np.pi / 180

def bearing(lat1, lon1, lat2, lon2):
    lat1 = lat1 * np.pi / 180
    lat2 = lat2 * np.pi / 180
    lon1 = lon1 * np.pi / 180
    lon2 = lon2 * np.pi / 180
    y = np.sin(lon2-lon1) * np.cos(lat2)
    x = np.cos(lat1)*np.sin(lat2) - np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)
    return np.arctan2(y, x)


# print(d13 * 180 / np.pi)

vlsra1, vlsrdec1 = (-5.6377748, 51.223348)
vlsra2, vlsrdec2 = (2.9147779, -23.337827)

# theta12 = bearing(run['ra'], run['dec'], run['ra']- run['pmra'], run['dec']- run['pmdec'])
theta12 = bearing(run['ra'].iloc[0], run['dec'].iloc[0], run['ra'].iloc[0]- vlsra1, run['dec'].iloc[0]- vlsrdec1)
theta13 = bearing(run['ra'].iloc[0], run['dec'].iloc[0], avgra, avgdec)
min_dist = np.arcsin(np.sin(d13) * np.sin(theta13-theta12))

print('MIN_DIST 1', min_dist[0] * 180 / np.pi)

# print(run['ra'], run['dec'], run['pmra'], run['pmdec'], run['parallax'])

# # pmlsr,ra, dec, pmra=pmra,pmdec=pmdec,par=par

#if vlsr conversion asks for radial velocities, set to some value

# c = SkyCoord(ra = run['ra'].iloc[:] * u.degree, dec = run['dec'].iloc[:]*u.degree, frame = 'icrs', pm_ra_cosdec = run['pmra'].iloc[:] *u.mas/u.yr, pm_dec = run['pmdec'].iloc[:] * u.mas/u.yr, distance = 1000/run['parallax'].iloc[:]*u.pc)
# l = 
# coords = ICRS(ra = run['ra'].iloc[0]*u.degree, dec = run['dec'].iloc[0]*u.degree, pm_ra_cosdec = run['pmra'].iloc[0]*u.mas/u.year, pm_dec = run['pmdec'].iloc[0]*u.mas/u.year)#, distance = 1/(run['parallax'].iloc[0] * u.mas))
# print(coords.proper_motion)
# lsr = coords.transform_to(LSR())
# print(lsr.proper_motion)

run['vlsrra'] = np.array([0,0])
run['vlsrdec'] = np.array([0,0])

c1 = SkyCoord(ra=run['ra'].iloc[0]*u.degree, dec=run['dec'].iloc[0]*u.degree, frame='icrs',pm_ra_cosdec=run['pmra'].iloc[0]*u.mas/u.yr, pm_dec=run['pmdec'].iloc[0]*u.mas/u.yr,distance=1000/run['parallax'].iloc[0]*u.pc,radial_velocity=0*u.km/u.s)
c2 = SkyCoord(ra=run['ra'].iloc[1]*u.degree, dec=run['dec'].iloc[1]*u.degree, frame='icrs',pm_ra_cosdec=run['pmra'].iloc[1]*u.mas/u.yr, pm_dec=run['pmdec'].iloc[1]*u.mas/u.yr,distance=1000/run['parallax'].iloc[1]*u.pc,radial_velocity=0*u.km/u.s)

l1=c1.transform_to(LSR())
# print(l1.pm_ra_cosdec.value)
# print(l1.pm_dec.value)
run['vlsrra'].iloc[0] = l1.pm_ra_cosdec.value
run['vlsrdec'].iloc[0] = l1.pm_dec.value

l2=c2.transform_to(LSR())
# print(l2.pm_ra_cosdec.value)
# print(l2.pm_dec.value)
run['vlsrra'].iloc[1] = l2.pm_ra_cosdec.value
run['vlsrdec'].iloc[1] = l2.pm_dec.value

theta12_vlsr = bearing(run['ra'].iloc[0], run['dec'].iloc[0], run['ra'].iloc[0]-run['vlsrra'].iloc[0], run['dec'].iloc[0]-run['vlsrdec'].iloc[0])
theta13_vlsr = bearing(run['ra'].iloc[0], run['dec'].iloc[0], avgra, avgdec)
min_dist_vlsr = np.arcsin(np.sin(d13) * np.sin(theta13_vlsr-theta12_vlsr))

print('MIN_DIST 2', min_dist_vlsr[0] * 180 / np.pi)

c_ori = SkyCoord(ra = avgra * u.degree, dec = avgdec * u.degree, frame = 'icrs', pm_ra_cosdec = np.mean(ori['pmra']) * u.mas/u.yr, pm_dec = np.mean(ori['pmdec'])*u.mas/u.yr, distance = 1000/np.mean(ori['parallax'])*u.pc, radial_velocity = 0 * u.km/u.s)
l_ori = c_ori.transform_to(LSR())
avgra_corr = avgra - (l_ori.pm_ra_cosdec.value / (3600 * 1000) * 2.6e6)
avgdec_corr = avgdec - (l_ori.pm_dec.value /(3600 * 1000) * 2.6e6)
print('avg vels', l_ori.pm_ra_cosdec.value, l_ori.pm_dec.value)

theta12_vlsr_ = bearing(run['ra'].iloc[0], run['dec'].iloc[0], run['ra'].iloc[0]-run['vlsrra'].iloc[0], run['dec'].iloc[0]-run['vlsrdec'].iloc[0])
theta13_vlsr_ = bearing(run['ra'].iloc[0], run['dec'].iloc[0], avgra_corr, avgdec_corr)
min_dist_vlsr_ = np.arcsin(np.sin(d13) * np.sin(theta13_vlsr_-theta12_vlsr_))

print('MIN_DIST 3', min_dist_vlsr_[0] * 180 / np.pi)

theta12_vlsr__ = bearing(run['ra'].iloc[0], run['dec'].iloc[0], run['ra'].iloc[0]-run['vlsrra'].iloc[0]-l_ori.pm_ra_cosdec.value, run['dec'].iloc[0]-run['vlsrdec'].iloc[0]-l_ori.pm_dec.value)
theta13_vlsr__ = bearing(run['ra'].iloc[0], run['dec'].iloc[0], avgra, avgdec)
min_dist_vlsr__ = np.arcsin(np.sin(d13) * np.sin(theta13_vlsr__-theta12_vlsr__))

print('MIN_DIST 4', min_dist_vlsr__[0] * 180 / np.pi)