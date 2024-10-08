from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from ejecutable import L
from astropy.io import ascii
from astropy.table import Table, hstack
import pandas as pd
import subprocess

Datos_L = L.group_by('index')
Galaxies = Datos_L.groups.keys

def crossmatch(galaxy):
    sex_catalog = ascii.read(f'sex/Extra/galaxy_{galaxy}')
    galaxies = Table.read(f'Catalog/Galaxies_Extra.csv')

    #Definición de las variables
    ra_sex = sex_catalog['ALPHA_SKY'].astype(float)
    dec_sex = sex_catalog['DELTA_SKY'].astype(float)

    ra_gal = galaxies['ra'].astype(float)
    dec_gal = galaxies['dec'].astype(float)

    catalog_sex = SkyCoord(ra=ra_sex, dec=dec_sex)
    catalog_gal = SkyCoord(ra=ra_gal*u.degree, dec=dec_gal*u.degree)

    #Crossmatch
    idx, d2d, d3d = match_coordinates_sky(catalog_sex, catalog_gal)
    
    #Radio de búsqueda
    max_sep = 3*u.arcsec
    sep_constraint = d2d < max_sep

    #Filtrar los datos
    matched1 = sex_catalog[sep_constraint]
    matched2 = galaxies[idx[sep_constraint]]

    combined_table = hstack([matched1, matched2])
    ra_com = combined_table['ra']
    if len(ra_com) > 1:
        print(f'Para la galaxia {galaxy}, se tomaron {len(ra_com)} fuentes.')
    else:
        print(f'Crossmatch correcto para la galaxia {galaxy}')
    combined_table.write(f'sex/Extra/Galaxy_{galaxy}.csv', format='csv', overwrite=True)
    return 

#ajustar = pd.read_csv('/home/seba/Documents/numeros_unicos.txt', header=None)
#n = ajustar[0].to_list()
for g in Galaxies['index']:
    crossmatch(g)
    



