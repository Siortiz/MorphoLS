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
from utils import filter_sel, coordenadas
import os
 
def median_sky(grupo, filtros):
    sk=[]
    for filtro in filtros:
        image = CCDData.read(f'/home/seba/Documents/DECALS/joined_bricks/{grupo}/{grupo}_image_{filtro}.fits', unit='adu')
        mean, median, std = sigma_clipped_stats(image.data, sigma=3.0)
        sk.append(median)
    texto_sky = ', '.join(str(sk))
    I = f'1) {texto_sky}'
    return I

def galfit_input(GRf, filtros, long):
    grupo=GRf['Group'][0]
    #Band labels
    f = filtros.tolist()
    texto_filtros = ', '.join(f)
    A1 = f"A1) {texto_filtros}"
    #Band wavelenghts
    w = long.tolist()
    texto_wvl = ', '.join(w)
    A2 = f"A2) {texto_wvl}"
    #Output data image block (FITS filename)
    B = 'B) galfitm_output/galfitm_group_{grupo}.fits'
    # Sigma image name (CSL of <nbands> FITS filenames or "none")
	# (if an individual filename is specified as "none", then that sigma
	#  image will be made from data; if the whole entry consists of just a
	#  single "none", then all sigma images will be made from data.)
	# One can also add a minimum sigma value, such that any galfit-created
	# sigma image will have a minimum of that value times the sky-subtracted
	# input data.
    C = 'C) none'
    # PSF fine sampling factor relative to data
    E = 'E) 1'
    #Mask
    mask = 'Field_Img/mask/mask_group_{grupo}.fits'
    num_filtros=len(filtros)
    F = 'F) '+', '.join([mask]*num_filtros)
    #File with parameter constraints (ASCII file)
    G = 'G) none'
    #Image region to fit (xmin xmax ymin ymax)
    H = 'H) 0, 1374, 0, 1374'
    #Size of the convolution box(x, y)
    I = 'I) 150, 150'
    #Magnitude photometric zeropoint
    J = 'J) '+', '.join(['22.5']*num_filtros)
    print(J)
    #Plate scale (dx dy) [arcsec per pixel]
    K = 'K) 0.262 0.262'
    #Display type (regular, curses, both)
    O = 'O) regular'
    #Options: 0=normal run; 1,2=make model/imgblock & quit
    P = 'P) 0'
    #Non-parametric component
    U = 'U) 0' #Standard parametric fitting
    W = 'W) input, model, residual'
    a = []
    d = []
    #Input data images (CSL of FITS filenames)
    for j in range(num_filtros):
        a.append(f'/home/seba/Documents/DECALS/joined_bricks/{grupo}/{grupo}_image_filtros[j].fits')
        d.append(f'Field_Img/psf/psf_group_{grupo}_filtros[j].fits')
    texto_img = ', '.join(a)
    A = 'A) {texto_img}'
    #Input PSF image (CSL of <nbands> FITS filenames
    texto_psf = ', '.join(d)
    D = 'D) {texto_psf}'
    Data = [A, A1, A2, B, C, D, E, F, G, H, I, J, K, O, P, U, W]
    for i in range(len(GRf)):
        Data.append(f'#----------Galaxia ({i})----------')
        Data.append('0) sersic') #Object type
        x = X[i]
        y = Y[i]
        Data.append('1) '+', '.join([str(x)]*num_filtros)+'1')
        Data.append('2) '+', '.join([str(y)]*num_filtros)+'1')
        #Magnitudes
        magnitudes = [GRf[f'mag_{filtro}'][i] for filtro in filtros]
        magnitudes_text = ', '.join(map(str, magnitudes))
        Data.append('3) '+magnitudes_text)
        
        #Filtrar sex_data
        mask = (sex_data['X_IMAGE'] == x) & (sex_data['Y_IMAGE'] == y)
        filtered_sex_data = sex_data[mask]

        R = filtered_sex_data['FLUX_RADIUS']
        Kron = filtered_sex_data['KRON_RADIUS']
        n = R/Kron
        texto_r = [str(R)]*num_filtros
        Data.append('4) '+', '.join(texto_r)+' 2')
        texto_n = [str(n)]
        Data.append('5) '+', '.join(texto_n)+' 2')
        el = 1/sex_data['ELONGATION']
        texto_el = [str(el)]*num_filtros
        Data.append('9) '+', '.join(texto_el)+' 1')
        if filtered_sex_data['THETA_IMAGE'] >= 0:
            Th = (filtered_sex_data['THETA_IMAGE']-90)
        else:
            Th = (filtered_sex_data['THETA_IMAGE']+90)
        texto_theta = [str(Th)]*num_filtros
        Data.append('10) '+', '.join(texto_theta)+' 1')
        Data.append('Z) 0') #Skip this model in output image? (yes=1, no=0)
        Data.append('#-------sky-------')
        Data.append('0) sky')
        Data.append(median_sky(grupo, filtros)) #sky background [ADU counts]
        Data.append('2) 0.000   0') #dsky/dx (sky gradient in x)	Data.append('2) 0.000      0 ') # dsky/dx (sky gradient in x) 
        Data.append('3) 0.000   0 ') # dsky/dy (sky gradient in y)
        Data.append('Z) 0')   # Skip this model in output image?  (yes=1, no=0)
        
        #Guardar cada línea de data en un archivo
        fic = open(f'inputs/galfit_{grupo}.input', 'w')
        for line in Data:
            print(line, file = fic)
        fic.close()
        return

L = Table.read('/home/seba/Documents/DECALS/Galaxies/Galaxies_DECALS_186.csv')
Datos_L = L.group_by('Group')
GL = Datos_L.groups.keys

for g in GL['Group']:
    filtros = filter_sel(g)[0]
    long = filter_sel(g)[1]
    X = coordenadas(g)[0]
    Y = coordenadas(g)[1]
    sex_data = ascii.read(f'sex/group_{g}')
    mask = Datos_L.groups.keys['Group'] == g
    GRf = Datos_L.groups[mask]
    galfit_input(GRf,filtros,long) 
