from astropy.table import *
from astropy.io import fits
import os
import numpy as np
L = Table.read('/home/seba/Documents/DECALS/Galaxies/Galaxies_DECALS_186.csv')
Datos_L = L.group_by('Group')
GL = Datos_L.groups.keys
        

def img_det(grupo):
    ruta_base = f'/home/seba/Documents/DECALS/joined_bricks/{grupo}'
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
        archivo = f'{ruta_base}/{grupo}_image_{filtro}.fits'
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
        im_det.writeto(f'Field_Img/det/det_group_{grupo}.fits', overwrite=True)
        print(f"Imagen de detección creada para el grupo {grupo} utilizando los filtros: {', '.join(filtros_usados)}.")
    else:
        print(f"No se encontraron archivos de los filtros especificados para crear la imagen de detección para el grupo {grupo}.")

    return

Data=[]
for g in GL['Group']:
    img_det(g)
    Data.append(f'sex Field_Img/det/det_group_{g}.fits -c sex.conf -CATALOG_NAME sex/group_{g} -CATALOG_TYPE ASCII_HEAD -PARAMETERS_NAME ./sex.param -DETECT_THRESH 3 -ANALYSIS_THRESH 3 -FILTER_NAME gauss_5.0_9x9.conv -SATUR_LEVEL 25000 -MAG_ZEROPOINT 22.5 -PIXEL_SCALE 0.262 -SEEING_FWHM 1.38 -CHECKIMAGE_TYPE SEGMENTATION -CHECKIMAGE_NAME Field_Img/det/det_group_{g}_seg.fits')


fic = open('sex_seg.sh', 'w')
for line in Data:
    print(line, file=fic)
fic.close
