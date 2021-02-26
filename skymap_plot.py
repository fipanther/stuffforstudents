import os
import numpy as np
import matplotlib.pyplot as plt

from ligo.skymap.io import fits
from ligo.skymap import postprocess

import healpy as hp

from astropy import units as u
from astropy.coordinates import SkyCoord


# this just makes a list of the files that need to be read to make the skymap
# I am reading through a folder that contains all the skymaps for GWTC-1
CWD = os.getcwd()
gwtc1 = CWD+'/gwtc1-skymaps'
#gwtc2 = CWD+'/all_skymaps'
all_files = os.listdir(gwtc1)
fnames = []
for file in all_files:
    if file.endswith('.fits'):
        fnames.append(file)



#can remove the loop if you just want one file, in that case you give the filename instead of the loop
# ras = []
# decs = []
for fname in fnames:
    print('processing '+fname+'...')
    fits_file = gwtc1+'/'+fname
    skymap, metadata = fits.read_sky_map(fits_file, nest=None)
    cls = 100 * postprocess.find_greedy_credible_levels(skymap)
    contour_levels = [50, 90, 99]

    #this was for Jamie trying to find the 'center' of the contours, ignore it
    # mx = np.max(skymap)
    # mean = np.mean(skymap)
    # std = np.std(skymap)
    # idx = np.where(skymap == mx)

    # print(mx, idx)
    # print('mean: {}'.format(mean))
    # print('std: {}'.format(std))

    # dec, ra= IndexToDeclRa(idx)
    # print('ra=',ra, 'dec=',dec)
    # ras.append(np.radians(ra[0]))
    # decs.append(np.radians(dec[0]))

    plt.figure(figsize=(15, 15))
    ax = plt.axes(projection='astro hours mollweide')
    cs = ax.contourf_hpx(
            (cls, 'ICRS'), nested=metadata['nest'], cmap='cylon_r')
    ax.grid()
    #if you want to plot a point somewhere on the projection, use this
    #ax.plot_coord(SkyCoord(ra,dec, unit='deg'), 'o',markeredgecolor='black', markersize=5)
    ax.tick_params(axis='both', which='major', labelsize=18)
    plt.savefig(fname+'.png')
    
    plt.show()