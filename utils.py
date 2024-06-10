import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy.table import Table
from ejecutable import L
#L = Table.read('/home/seba/Documents/DECALS/Galaxies/Galaxies_DECALS_186.csv')
Datos_L = L.group_by('Group')
GL = Datos_L.groups.keys


def filter_sel(grupo):
    filtros_existentes=[]
    lambda_por_filtro=[]
    filtros_longitudes = {'g': '4770', 'r': '6231', 'i': '7625', 'z': '9134'}
    
    directorio = f'/home/seba/Documents/DECALS/joined_bricks/{grupo}/'
    if os.path.exists(directorio):
        for filtro, lambda_ in filtros_longitudes.items():
            nombre_archivo = f'{grupo}_image_{filtro}.fits'
            if  nombre_archivo in os.listdir(directorio):
                filtros_existentes.append(filtro)
                lambda_por_filtro.append(lambda_)
        filtros_existentes = np.array(filtros_existentes)
        lambda_por_filtro = np.array(lambda_por_filtro)
        return filtros_existentes, lambda_por_filtro
    else:
        print(f'El grupo {grupo} no existe')
def coordenadas(grupo):
    #Coordenadas a posiciones
    hdu = fits.open(f'Field_Img/det/det_group_{grupo}.fits')
    hdur = hdu[0].header
    wcs = WCS(header=hdu[0].header)
    X_0 = hdur['CRPIX1']
    Y_0 = hdur['CRPIX2']
    RA_0 = hdur['CRVAL1']*u.degree
    DEC_0 = hdur['CRVAL2']*u.degree
    delta = np.abs(hdur['CDELT1'])*u.degree
    X = []
    Y = []
    for i in range(len(Datos_L)):
        RA = Datos_L['ra'][i]*u.degree
        DEC = Datos_L['dec'][i]*u.degree
        (X_new, Y_new) = wcs.all_world2pix(RA, DEC, 1)
        X.append(X_new)
        Y.append(Y_new)
    X=np.array(X)
    Y=np.array(Y)
    print(X, Y)
    return X, Y
