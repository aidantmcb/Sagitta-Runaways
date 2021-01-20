import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import astropy.units as u
from tol_colors import tol_cmap

PRGn = tol_cmap("sunset")

from astropy.coordinates import SkyCoord

def plotSky(tab, aux, avg_l, avg_b):
    fig = plt.figure(figsize = (12,9))
    ax = plt.subplot(projection = 'aitoff')
    ax.grid(True)

    coords = SkyCoord(-tab['l'], tab['b'], unit = 'deg', frame='galactic')
    c_l = coords.l.wrap_at(180*u.deg).radian
    c_b = coords.b.radian

    avgs = SkyCoord(-avg_l, avg_b, unit = 'deg', frame = 'galactic')
    avg_l_rad = avgs.l.wrap_at(180*u.deg).radian
    avg_b_rad = avgs.b.radian
    a = ax.scatter(c_l, c_b, s=0.05, c = aux, cmap = PRGn)
    
    b = ax.scatter(avg_l_rad, avg_b_rad, c = 'orange', marker = "*")

    fig.colorbar(a, ax = ax, label = 'Aux')
    plt.show()