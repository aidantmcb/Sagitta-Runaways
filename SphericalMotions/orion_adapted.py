import numpy as np
from numpy import arctan2 as atan2, sin,cos,arccos as acos,sqrt
from astropy import units as u
from astropy.coordinates import SkyCoord, LSR

ae_c = SkyCoord(ra=79.0756206*u.degree, dec=34.31231862*u.degree, frame='icrs',pm_ra_cosdec=-4.440*u.mas/u.yr, pm_dec=43.368*u.mas/u.yr,distance=1000/2.4642*u.pc,radial_velocity=0*u.km/u.s)
ae_l=ae_c.transform_to(LSR())
mu_c = SkyCoord(ra=86.49956401*u.degree, dec=-32.30643465*u.degree, frame='icrs',pm_ra_cosdec=2.988*u.mas/u.yr, pm_dec=-22.030*u.mas/u.yr,distance=1000/2.1476*u.pc,radial_velocity=0*u.km/u.s)
mu_l=mu_c.transform_to(LSR())

cluster = SkyCoord(ra=83.8*u.degree, dec=-5.4*u.degree, frame='icrs',pm_ra_cosdec=0.8*u.mas/u.yr, pm_dec=2.8*u.mas/u.yr,distance=390*u.pc,radial_velocity=0*u.km/u.s)
star=ae_l

lon1=star.ra.value*np.pi/180
lat1=star.dec.value*np.pi/180
lon2=cluster.ra.value*np.pi/180
lat2=cluster.dec.value*np.pi/180
lon3=(star.ra.value-(star.pm_ra_cosdec.value-cluster.pm_ra_cosdec.value)/cos(star.dec.value*np.pi/180)/1e3/3600)*np.pi/180
lat3=(star.dec.value-(star.pm_dec.value-cluster.pm_dec.value)/1e3/3600)*np.pi/180

bearing1=atan2(cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(lon2-lon1),sin(lon2-lon1)*cos(lat2))
bearing2=atan2(cos(lat1)*sin(lat3)-sin(lat1)*cos(lat3)*cos(lon3-lon1),sin(lon3-lon1)*cos(lat3))
dist=acos( sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(lon2-lon1) )/np.pi*180
velocity=sqrt((star.pm_ra_cosdec.value-cluster.pm_ra_cosdec.value)**2+(star.pm_dec.value-cluster.pm_dec.value)**2)


print('Angle with proper motions being off by:', (bearing1-bearing2)/np.pi*180, 'deg')
print('Radial alignment:', cos((bearing1-bearing2)))
print('Distance:',dist,'deg')
print('Velocity:',velocity,'mas/yr')
print('Traceback time:', dist/(velocity/1e3/3600*1e6),'Myr')