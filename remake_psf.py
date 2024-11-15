from astropy.io import ascii, fits
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from astropy.nddata import NDData, CCDData
from photutils.psf import extract_stars, EPSFBuilder
from astropy.wcs import WCS
from utils import filter_sel
from ejecutable import L
import numpy as np
import astropy.units as u
import os



#Desenmascarar fuentes
det = f'Field_Img/det/det_group_35_seg.fits'
segment_labels = [94, 185, 246, 277]
fig = fits.open(det)
    
for label in segment_labels:
    arr_mask = fig[0].data == label     
    fig[0].data[arr_mask] = 0
    
arr_mask = fig[0].data > 0
fig[0].data[arr_mask] = 1
hdr = fig[0].header
imgF = fits.PrimaryHDU(fig[0].data, header=hdr)
imgF.writeto(f'Field_Img/mask/mask_group_35.fits', overwrite=True)

