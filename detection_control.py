from astropy.table import *
from astropy.io import fits
import os
import numpy as np
from utils import filter_sel 
from ejecutable import L
#L = Table.read('Catalog/GalfitM_DECALS_control_sample.csv')
Datos_L = L.group_by('index')
GL = Datos_L.groups.keys
        

def img_det(galaxy):
    ruta_base = f'/home/seba/Documents/DECALS/joined_bricks_cs/{galaxy}'
    filtros_prioritarios = ['g', 'r', 'z']
    filtros_secundarios = ['i']
    suma_imagen = None
    hdr = None

    # Lista de todos los filtros en orden de prioridad
    todos_filtros = filtros_prioritarios + filtros_secundarios
    
    #Lista para almacenar los filtros utilizados
    filtros_usados = []

    # Iterar sobre los filtros en orden de prioridad
    for filtro in todos_filtros:
        archivo = f'{ruta_base}/{galaxy}_image_{filtro}.fits'
        if os.path.exists(archivo):
            with fits.open(archivo) as hdulist:
                imagen = hdulist[0].data
                if suma_imagen is None:
                    suma_imagen = np.zeros_like(imagen)
                    hdr = hdulist[0].header  # Usar el header del primer archivo encontrado
                suma_imagen += imagen
                filtros_usados.append(filtro)
            # Si ya hemos agregado los tres filtros prioritarios, no necesitamos más
            if len(filtros_usados) == 3 and all(f in filtros_usados for f in filtros_prioritarios):
                break

    if suma_imagen is not None:
        im_det = fits.PrimaryHDU(suma_imagen, header=hdr)
        im_det.writeto(f'Field_Img/det_cs/det_galaxy_{galaxy}.fits', overwrite=True)
        print(f"Imagen de detección creada para la galaxia {galaxy} utilizando los filtros: {', '.join(filtros_usados)}.")
    else:
        print(f"No se encontraron archivos de los filtros especificados para crear la imagen de detección para la galaxia {galaxy}.")

    return

def fwhm(GR, filtros):
    psfsize = []
    for filtro in filtros:
        mean_band = np.mean(GR[f'psfsize_{filtro}'])
        #print(GR[f'psfsize_{filtro}'])
        psfsize.append(mean_band)
    mean_fwhm = np.mean(psfsize)
    return mean_fwhm
Data=[]
for g in GL['index']:
    if g > 1040:
        try:
            img_det(g)
            print(f'Galaxia {g}')
        except:
            print(f'No se pudo crear la imagen de detección de la galaxia {g}')
        filtros = filter_sel(g)[0]
        mask = Datos_L.groups.keys['index']==g
        GR = Datos_L.groups[mask]
        #print(filtros)
        #print(GR['psfsize_g'], GR['psfsize_i'], GR['psfsize_z'])
        fwhm_size = fwhm(GR, filtros)
        Data.append(f'sex Field_Img/det_cs/det_galaxy_{g}.fits -c sex.conf -CATALOG_NAME sex/galaxy_{g} -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME ./sex.param -DETECT_THRESH 1 -ANALYSIS_THRESH 1 -FILTER_NAME gauss_5.0_9x9.conv -SATUR_LEVEL 25000 -MAG_ZEROPOINT 22.5 -PIXEL_SCALE 0.262 -SEEING_FWHM {fwhm_size} -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME Field_Img/det_cs/det_galaxy_{g}_seg.fits')


fic = open('sex_seg.sh', 'w')
for line in Data:
    print(line, file=fic)
fic.close
