import numpy as np
import pandas as pd
from astropy.io import fits
import os

fig_directory = 'c:/users/sahal/Desktop/Plots/b/Cluster-'
fname = 'c:/Users/sahal/Desktop/RunawayDetector_2-10_Taurus.fits'
data = np.array(fits.open(fname)[1].data)
df = pd.DataFrame(data.byteswap().newbyteorder())

labels = np.sort(np.unique(df['cluster_label']))

def get_dist(l, b, avg_l, avg_b):
        d = np.sqrt((l-avg_l)**2 + (b-avg_b)**2)
        return d

for label in [11]:
    label = int(label)
    print(label)
    c = df.iloc[np.where(df['cluster_label'] == label)[0]]
    avg_l, avg_b = np.unique(c['avg_l'])[0], np.unique(c['avg_b'])[0]
    max_dist = np.max(get_dist(c['l'], c['b'], c['avg_l'], c['avg_b']))
    max_pmoff = np.max(get_dist(c['vlsrl'], c['vlsrb'], c['avg_pml'], c['avg_pmb'] ))
    radius = np.min([max_dist, 45]) 
    title = 'Label ' + str(label) + ' ID ' + str(np.unique(c['cluster_id']))

    fig_path = fig_directory + str(label) + '.jpg'
    STILTS = r'''java -jar topcat-lite.jar -stilts plot2sky xpix=600 ypix=600 sex=false title="{title}" fontsize=20 clon={avg_l} clat={avg_b} radius={r} auxmap=plasma auxvisible=true auxlabel=Traceback legend=false in={input} icmd="select 'cluster_label == {label}'" layer_1=Mark lon_1=l lat_1=b aux_1="exp10(traceback - 6) " shading_1=aux layer_2=SkyVector lon_2=l lat_2=b dlon_2=vlsrl-avg_pml dlat_2=vlsrb-avg_pmb aux_2="exp10(traceback-6)" shading_2=aux scale_2=0.09 layer_3=Mark lon_3=avg_l lat_3=avg_b shading_3=auto shape_3=filled_square size_3=4 out={outpath} ofmt=jpeg forcebitmap=true'''.format(
        title = title, avg_l = avg_l, avg_b = avg_b, r = radius, input = fname, label = label, outpath = fig_path
)
    os.chdir('C:/Users/Sahal/Desktop')
    os.system(STILTS)