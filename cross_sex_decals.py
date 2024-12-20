from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky
from ejecutable import L
from astropy.io import ascii
from astropy.table import Table, hstack
import pandas as pd
import subprocess

Datos_L = L.group_by('Group')
Grupos = Datos_L.groups.keys

def crossmatch(grupo):
    sex_catalog = ascii.read(f'sex/Groups/group_{grupo}')
    galaxies = Table.read(f'/home/seba/Documents/DECALS/Galaxies/Galaxies_DECALS_{grupo}.csv')

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
    max_sep = 5*u.arcsec
    sep_constraint = d2d < max_sep

    #Filtrar los datos
    matched1 = sex_catalog[sep_constraint]
    matched2 = galaxies[idx[sep_constraint]]

    combined_table = hstack([matched1, matched2])
    ra_com = combined_table['ra']
    print(matched1)
    print(matched2)
    print(galaxies)
    print(f'La cantidad de galaxias en el grupo {grupo} antes del crossmatch es {len(ra_gal)} y después del crossmatch es {len(ra_com)}')
    combined_table.write(f'sex/Groups/Galaxies_group_{grupo}.csv', format='csv', overwrite=True)
    return 

crossmatch(193)
#for g in Grupos['Group']:
#    if g != 275:
#        crossmatch(g)
    



