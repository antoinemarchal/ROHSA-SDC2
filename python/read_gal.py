import numpy as np
from astropy.io import fits   
import matplotlib.pyplot as plt

path = "../output/"
fitsname = "sky_ldev_v2_1_data.fits"
hdu = fits.open(path+fitsname) 
hdr = hdu[0].header 
cube = hdu[0].data
