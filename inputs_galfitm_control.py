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
from ejecutable import L
import os
import math
import pandas as pd
def median_sky(galaxy, filtros):
    sk=[]
    for filtro in filtros:
        image = CCDData.read(f'/home/seba/Documents/DECALS/joined_bricks_cs/{galaxy}/{galaxy}_image_{filtro}.fits', unit='adu')
        mean, median, std = sigma_clipped_stats(image.data, sigma=3.0)
        sk.append(median)
    #print(sk)
    texto_sky = ','.join(map(str,sk))
    I = f'1) {texto_sky}'
    return I

def galfit_input(GRf, filtros, long, size):
    galaxy=GRf['index'][0]
    #Band labels
    f = filtros.tolist()
    texto_filtros = ','.join(f)
    A1 = f"A1) {texto_filtros}"
    #Band wavelenghts
    w = long.tolist()
    texto_wvl = ','.join(w)
    A2 = f"A2) {texto_wvl}"
    #Output data image block (FITS filename)
    B = f'B) galfitm_output/Control_Sample/galfitm_galaxy_{galaxy}.fits'
    # PSF fine sampling factor relative to data
    E = 'E) 1'
    #Mask
    mask = f'Field_Img/mask_cs/mask_galaxy_{galaxy}.fits'
    num_filtros=len(filtros)
    F = 'F) '+','.join([mask]*num_filtros)
    #File with parameter constraints (ASCII file)
    G = 'G) none'
    #Image region to fit
    H = f'H) 0 {size} 0 {size}'
    #Size of the convolution box(x, y)
    I = 'I) 300 300'
    #Magnitude photometric zeropoint
    J = 'J) '+','.join(['22.5']*num_filtros)
    #Plate scale (dx dy) [arcsec per pixel]
    K = 'K) 0.262 0.262'
    #Display type (regular, curses, both)
    O = 'O) regular'
    #Options: 0=normal run; 1,2=make model/imgblock & quit
    P = 'P) 0'
    #Non-parametric component
    U = 'U) 0' #Standard parametric fitting
    V = 'V) 0' #Multinest
    W = 'W) input,model,residual'
    a = []
    c = []
    d = []
    #Input data images (CSL of FITS filenames)
    for j in range(num_filtros):
        a.append(f'/home/seba/Documents/DECALS/joined_bricks_cs/{galaxy}/{galaxy}_image_{filtros[j]}.fits')
        c.append(f'/home/seba/Documents/DECALS/sigma_image_cs/{galaxy}/sigma_{galaxy}_weight_{filtros[j]}.fits')
        d.append(f'Field_Img/psf_cs/psf_galaxy_{galaxy}_{filtros[j]}.fits')       
    texto_img = ','.join(a)
    A = f'A) {texto_img}'
    # Sigma image name (CSL of <nbands> FITS filenames or "none")
    # (if an individual filename is specified as "none", then that sigma
    #  image will be made from data; if the whole entry consists of just a
    #  single "none", then all sigma images will be made from data.)
    # One can also add a minimum sigma value, such that any galfit-created
    # sigma image will have a minimum of that value times the sky-subtracted
    # input data.
    texto_sigma = ','.join(c)
    C = f'C) {texto_sigma}'    
    #Input PSF image (CSL of <nbands> FITS filenames
    texto_psf = ','.join(d)
    D = f'D) {texto_psf}'
    Data = [A, A1, A2, B, C, D, E, F, G, H, I, J, K, O, P, U, W]
    X_1=[]
    Y_1=[]
    for i in range(len(sex_data)):
        Data.append(f'#----------Galaxia ({i})----------')
        Data.append('0) sersic') #Object type
        x = X[i]
        y = Y[i]
        X_1.append(x)
        Y_1.append(y)
        Data.append('1) '+','.join([str(x)]*num_filtros)+' 1')
        Data.append('2) '+','.join([str(y)]*num_filtros)+' 1')
        #Magnitudes
        magnitudes = [GRf[f'mag_{filtro}'][i] for filtro in filtros]
        # Reemplaza valores infinitos por 20
        magnitudes = [20 if math.isinf(mag) else mag for mag in magnitudes]
        magnitudes_text = ','.join(map(str, magnitudes))
        Data.append('3) '+magnitudes_text+f' {num_filtros}')
        
                
        R = sex_data['FLUX_RADIUS'][i]
        Kron = sex_data['KRON_RADIUS'][i]
        n = R/Kron
        texto_r = [str(R)]*num_filtros
        Data.append('4) '+','.join(texto_r)+' 2')
        texto_n = [str(n)]*num_filtros
        Data.append('5) '+','.join(texto_n)+' 2')
        el = 1/sex_data['ELONGATION'][i]
        texto_el = [str(el)]*num_filtros
        Data.append('9) '+','.join(texto_el)+' 1')
        if sex_data['THETA_IMAGE'][i] >= 0:
            Th = (sex_data['THETA_IMAGE'][i] - 90)
        else:
            Th = (sex_data['THETA_IMAGE'][i] + 90)
        texto_theta = [str(Th)]*num_filtros
        Data.append('10) '+','.join(texto_theta)+' 1')
    Data.append('Z) 0') #Skip this model in output image? (yes=1, no=0)
    Data.append('#-------sky-------')
    Data.append('0) sky')
    Data.append(median_sky(galaxy, filtros)+' 1') #sky background [ADU counts]
    Data.append('2) 0.000   0') #dsky/dx (sky gradient in x)	Data.append('2) 0.000      0 ') # dsky/dx (sky gradient in x) 
    Data.append('3) 0.000   0 ') # dsky/dy (sky gradient in y)
    Data.append('Z) 0')   # Skip this model in output image?  (yes=1, no=0)
        
    #Guardar cada línea de data en un archivo
    fic = open(f'inputs/Control_Sample/galfit_{galaxy}.input', 'w')
    for line in Data:
        print(line, file = fic)
    fic.close()

    print(f'Archivo de input generado para la galaxia {galaxy}')
    return X_1, Y_1

#L = Table.read('/home/seba/Documents/DECALS/Galaxies/Galaxies_DECALS_186.csv')
Datos_L = L.group_by('index')
GL = Datos_L.groups.keys
ajustar = pd.read_csv('/home/seba/Documents/numeros_unicos.txt', header=None)
n = ajustar[0].to_list()
for g in GL['index']:
    if g in n and g > 204:
        with fits.open(f'/home/seba/Documents/DECALS/joined_bricks_cs/{g}/{g}_image_g.fits') as hdul:
            size = hdul[0].data.shape[0]
        filtros = filter_sel(g)[0]
        long = filter_sel(g)[1]
        sex_data = Table.read(f'sex/Control_Sample/Galaxy_{g}.csv')
        X = sex_data['X_IMAGE']
        Y = sex_data['Y_IMAGE']
        mask = Datos_L.groups.keys['index'] == g
        GRf = Datos_L.groups[mask]
        X_1, Y_1 = galfit_input(GRf,filtros,long, size) 
