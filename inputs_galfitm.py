#Genera el archvo de input para galfitm
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import *
from astropy.table import Table
from astropy.io import ascii
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
from astropy.nddata import CCDData
from astropy.table import *
from utils import filter_sel

 
def median_sky(grupo):
    sk=[]
    for filtro in filtros:
        image = CCDData.read(f'/home/seba/Documents/DECALS/joined_bricks/{grupo}/{grupo}_image_{filtro}', unit='adu')
        if os.path.exists(image):
            mean, median, std = sigma_clipped_stats(image.data, sigma=3.0)
            sk.append(median)
        return sk

def galfit_input(grupo, filtros, long):
    #Band labels
    f = filtros.tolist()
    texto_filtros = ', '.join(f)
    A1 = f"A1) {texto_filtros}"
    #Band wavelenghts
    w = long.tolist()
    texto_wvl = ', '.join(w)
    A2 = f"A2) {texto_wvl}"
    #Output data image block (FITS filename)
    B = 'B) galfitm_output/galfitm_group_{group}.fits'
    # Sigma image name (CSL of <nbands> FITS filenames or "none")
	# (if an individual filename is specified as "none", then that sigma
	#  image will be made from data; if the whole entry consists of just a
	#  single "none", then all sigma images will be made from data.)
	# One can also add a minimum sigma value, such that any galfit-created
	# sigma image will have a minimum of that value times the sky-subtracted
	# input data.
	C='C) none'
    # PSF fine sampling factor relative to data
    E = 'E) 1'
    #Mask

L = Table.read('/home/seba/Documents/DECALS/Galaxies/Galaxies_DECALS_186.csv')
Datos_L = L.group_by('Group')
GL = Datos_L.groups.keys

for g in GL['Group']:
    filtros = filter_sel(g)[0]
    long = filter_sel(g)[1]
    print(filtros, long)
    galfit_input(g,filtros,long) 
