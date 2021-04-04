import numpy as np
from astropy.io import fits   
import matplotlib.pyplot as plt

path = "/gpfsdswork/dataset/SKA-DC2/"
fitsname = "sky_full_v2.fits"
hdu = fits.open(path+fitsname) 
hdr = hdu[0].header 
# cube = hdu[0].data

#Export header
empty_primary = fits.PrimaryHDU(header=hdr)
hdulist = fits.HDUList([empty_primary])
hdulist.writeto("sky_full_v2_hdr.fits", clobber=True)
