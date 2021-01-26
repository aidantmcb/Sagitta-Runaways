import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from tol_colors import tol_cmap

import subprocess
import sys
import os

PRGn = tol_cmap("sunset")

from astropy.coordinates import SkyCoord

def plotSky(tab, aux, avg_l, avg_b):
    fig = plt.figure(figsize = (12,9))
    ax = plt.subplot(projection = 'aitoff')
    ax.grid(True)

    coords = SkyCoord(-tab['l'], tab['b'], unit = 'deg', frame='galactic')
    c_l = coords.l.wrap_at(180*u.deg).radian
    c_b = coords.b.radian

    avgs = SkyCoord(avg_l, avg_b, unit = 'deg', frame = 'galactic')
    avg_l_rad = avgs.l.wrap_at(180*u.deg).radian
    avg_b_rad = avgs.b.radian
    a = ax.scatter(-c_l, c_b, s=0.25, c = aux, cmap = PRGn)
    
    b = ax.scatter(-avg_l_rad, avg_b_rad, c = 'orange', marker = "*")

    fig.colorbar(a, ax = ax, label = 'Aux')
    plt.show()

def plotSTILTS(avg_l, avg_b, avg_pml, avg_pmb, min_yso, max_age, radius = 35, align_threshold = 0.2, l = 'l'):
    align = r'''cos(atan2(vlsrl-{avg_pml},vlsrb-{avg_pmb})-atan2({lat}-{avg_l},b-{avg_b}))'''.format(
        avg_pml = avg_pml, avg_pmb = avg_pmb, lat = l, avg_l = avg_l, avg_b = avg_b
    )
    pm_mag = r'''sqrt(pow(vlsrl-{avg_pml},2)+pow(vlsrb-{avg_pmb},2))'''.format(
        avg_pml = avg_pml, avg_pmb = avg_pmb
    ) #in mas/yr
    velo = pm_mag + '*1000/parallax*(3.086e13)/(3.154e7)/13600000*3.14/180' # in km/s
    traceback = r'''sqrt(pow({lat}-{avg_l},2)+pow(b-{avg_b},2))/'''.format(
        lat = l, avg_l = avg_l, avg_b = avg_b, 
    ) + velo
    STILTS_str = r'''java -jar topcat-lite.jar -stilts plot2sky xpix=1161 ypix=869 scalebar=false sex=false texttype=latex fontsize=24 fontstyle=serif clon={clon} clat={clat} radius={radius} auxmap=Pastel auxflip=false auxvisible=true auxlabel="Aux" legend=false in="C:\Users\sahal\Desktop\sagitta_edr3-l1.fits" icmd="select 'yso > {yso} & age < {age} & {alignment} > 1-{align_threshold} & {velocity} > {velocity_min}'" lon=L lat=B aux={aux} shading=aux layer_1=Mark size_1=2 layer_2=SkyVector dlon_2=vlsrl-{avg_pml} dlat_2=vlsrb-{avg_pmb} scale_2=0.07'''.format(
        clon = avg_l, clat = avg_b, radius = radius, yso = min_yso, age = max_age, alignment = align, align_threshold = align_threshold, 
        velocity = velo, velocity_min = 10, avg_pml = avg_pml, avg_pmb = avg_pmb, aux = traceback)

    os.chdir('C:/Users/Sahal/Desktop')
    os.system(STILTS_str)
    print(STILTS_str)


