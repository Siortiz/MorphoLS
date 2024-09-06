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
from tqdm import tqdm
import pandas as pd
#L = Table.read('/home/seba/Documents/DECALS/Galaxies/Galaxies_DECALS_5.csv')
Datos_L = L.group_by('index')
GL = Datos_L.groups.keys
#filtros = filter_sel(Datos_L['Group'][0])
#print(filtros)
def psf_maker(galaxy, filtro):
    hdu = fits.open(f'/home/seba/Documents/DECALS/joined_bricks_cs/{galaxy}/{galaxy}_image_{filtro}.fits')
    data=hdu[0].data
    source_catalog = ascii.read(f'/home/seba/Documents/MorphoLS/sex/galaxy_{galaxy}')
    stars = source_catalog[(source_catalog['FLUX_RADIUS'] < 3) & (source_catalog['MAG_AUTO'] < 22)]
    x_sex = source_catalog['X_IMAGE']
    y_sex = source_catalog['Y_IMAGE']
    x = stars['X_IMAGE']
    y = stars['Y_IMAGE']
    snr = stars['SNR_WIN']
    #Para no considerar las estrellas de los bordes
    
    size_box = 25
    hsize = (size_box) - 1/2
    mask = ((x > hsize) & (x < (data.shape[1] -1 - hsize)) & (y > hsize) & (y <   (data.shape[0] -1 - hsize)))
    x=x[mask]
    y=y[mask]
    snr=snr[mask]
    selected_x = []
    selected_y = []
    selected_snr = []
    # Itera sobre cada estrella y verifica si hay otras estrellas dentro de la misma caja
    for i in range(len(x)):
    # Coordenadas de la estrella actual
        x_star = x[i]
        y_star = y[i]
        snr_star = snr[i]
    
    # Define los límites de la caja alrededor de la estrella
        x_min = int(x_star - hsize)
        x_max = int(x_star + hsize)
        y_min = int(y_star - hsize)
        y_max = int(y_star + hsize)
    
    # Verifica si hay otras estrellas dentro de la caja
        other_stars_inside = any((x_sex[j] >= x_min and x_sex[j] <= x_max and y_sex[j] >= y_min and y_sex[j] <= y_max) for j in range(len(x_sex)) if x_sex[j] != x[i])
    
    # Si no hay otras estrellas dentro de la caja, agrega la estrella actual a las seleccionadas
        if not other_stars_inside:
            selected_x.append(x_star)
            selected_y.append(y_star)
            selected_snr.append(snr_star)

    # Crea una nueva tabla con las estrellas seleccionadas
    stars_tbl = Table()
    stars_tbl['x'] = selected_x
    stars_tbl['y'] = selected_y
    stars_tbl['SNR'] = selected_snr

    #print(f'La psf se realizará con {len(stars_tbl)} estrellas')
    stars_tbl = stars_tbl[(stars_tbl['SNR'] > 750) & (stars_tbl['SNR'] < 25000)]
    print(f'La psf se realizará con {len(stars_tbl)} estrellas')
    #Extraer el background de la imagen
    
    image = CCDData.read(f'/home/seba/Documents/DECALS/joined_bricks_cs/{galaxy}/{galaxy}_image_{filtro}.fits', unit='adu')
    mean_val, median_val, std_val = sigma_clipped_stats(image.data, sigma=3.0)
    data -= median_val
    nddata = NDData(data=data)

    #Extraer estrellas
    
    extracted_stars = extract_stars(nddata, stars_tbl, size=27)
    
    #Calcular la psf
    
    epsf_builder = EPSFBuilder(oversampling=1, maxiters=3, progress_bar=True)
    epsf, fitted_stars = epsf_builder(extracted_stars)

    hdu_psf = fits.PrimaryHDU(epsf.data)
    hdul_psf = fits.HDUList([hdu_psf])
    hdul_psf.writeto(f'Field_Img/psf_cs/psf_galaxy_{galaxy}_{filtro}.fits', overwrite=True)
    
    return

def mask(galaxy):
    #Coordenadas de las galaxias

    Gal_se = Table.read(f'sex/Control_Sample/Galaxy_{galaxy}.csv')

    X=np.array(Gal_se['X_IMAGE'])
    Y=np.array(Gal_se['Y_IMAGE'])

    #Desenmascarar fuentes
    det = f'Field_Img/det_cs/det_galaxy_{galaxy}_seg.fits'
    if os.path.exists(det):
        fig = fits.open(det)
    
        for i in range(len(X)):
            arr_mask = fig[0].data == fig[0].data[int(Y[i])][int(X[i])]
         
            fig[0].data[arr_mask] = 0
        arr_mask = fig[0].data > 0
        fig[0].data[arr_mask] = 1
        hdr = fig[0].header
        imgF = fits.PrimaryHDU(fig[0].data, header=hdr)
        imgF.writeto(f'Field_Img/mask_cs/mask_galaxy_{galaxy}.fits', overwrite=True)
        return
    else:
        print('La imagen de detección de la galaxia {galaxy} no existe')
galaxias_no_psf = []
#new_try = pd.read_csv('/home/seba/Documents/numeros_unicos.txt', header=None)
#no_ajustar = pd.read_csv('/home/seba/Documents/numeros_unicos.txt', header=None)
#n = no_ajustar[0].to_list()
#print(n)
#for a in new_try[:0]:
    #print(a)
#print(new_try_list)
for g in GL['index']:
    filtros = filter_sel(g)[0]
    print(g)
    mask(g)
    for filtro in filtros:
        try:
            psf_maker(g, filtro)
        except:
            galaxias_no_psf.append(g)    
            print(f'La psf de la galaxia {g} no pudo ser calculada')
print(galaxias_no_psf)









