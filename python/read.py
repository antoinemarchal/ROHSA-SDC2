import numpy as np
from astropy.io import fits   
import matplotlib.pyplot as plt

plt.ion()

fitsname = "/mnt/raid-cita/amarchal/SDC2_10G/data/development/output/sky_dev_cat_1.fits" 
hdu = fits.open(fitsname) 
hdr = hdu[0].header 
cube = hdu[0].data 
