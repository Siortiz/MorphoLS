from astropy.table import *
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import os
import re
from astropy.table import Table
from utils import filter_sel
#from ejecutable import L

#fil_name = np.array(['g', 'r', 'i', 'z'])

def leer_header(NGAL, HEADER, fil_name, n_filtros):
    data = {}
    data['CHI2NU'] = HEADER.get('CHI2NU', 0)
    
    for filt in all_filters:
        if filt in fil_name:
            xc = HEADER.get(f'{NGAL}_XC_{filt}', '0.0 +/- 0.0').split()
            yc = HEADER.get(f'{NGAL}_YC_{filt}', '0.0 +/- 0.0').split()
            R = HEADER.get(f'{NGAL}_RE_{filt}', '0.0 +/- 0.0').split()
            m = HEADER.get(f'{NGAL}_MAG_{filt}', '0.0 +/- 0.0').split()
            N = HEADER.get(f'{NGAL}_n_{filt}', '0.0 +/- 0.0').split()
            ar = HEADER.get(f'{NGAL}_AR_{filt}', '0.0 +/- 0.0').split()
            PA = HEADER.get(f'{NGAL}_PA_{filt}', '0.0 +/- 0.0').split()
        else:
            xc = yc = R = m = N = ar = PA = ['0.0', '+/-', '0.0']

        # Comprobar la longitud de la lista antes de acceder a índices
        data[f'XC_{filt}'] = xc[0] if len(xc) > 0 else '0.0'
        data[f'e_XC_{filt}'] = xc[2] if len(xc) > 2 else '0.0'
        data[f'YC_{filt}'] = yc[0] if len(yc) > 0 else '0.0'
        data[f'e_YC_{filt}'] = yc[2] if len(yc) > 2 else '0.0'
        data[f'RE_{filt}'] = R[0] if len(R) > 0 else '0.0'
        data[f'e_RE_{filt}'] = R[2] if len(R) > 2 else '0.0'
        data[f'MAG_{filt}'] = m[0] if len(m) > 0 else '0.0'
        data[f'e_MAG_{filt}'] = m[2] if len(m) > 2 else '0.0'
        data[f'n_{filt}'] = N[0] if len(N) > 0 else '0.0'
        data[f'e_n_{filt}'] = N[2] if len(N) > 2 else '0.0'
        data[f'AR_{filt}'] = ar[0] if len(ar) > 0 else '0.0'
        data[f'e_AR_{filt}'] = ar[2] if len(ar) > 2 else '0.0'
        data[f'PA_{filt}'] = PA[0] if len(PA) > 0 else '0.0'
        data[f'e_PA_{filt}'] = PA[2] if len(PA) > 2 else '0.0'

    return data

def grafico(fits, galaxy, n_filtros):
    fig, axs = plt.subplots(3, 4, figsize=(20, 10), sharex=True, sharey=True)
    for j in range(n_filtros):
        axs[0, j].imshow(fits[j].data, cmap='gray', vmin=-0.5, vmax=0.5)
        axs[0, j].set_title(f'{fil_name[j]}')
        axs[1, j].imshow(fits[n_filtros + j].data, cmap='gray', vmin=-0.5, vmax=0.5)
        axs[2, j].imshow(fits[2 * n_filtros + j].data, cmap='gray', vmin=-0.5, vmax=0.5)
        axs[0, 0].text(60, 140, '3 arcmin', fontsize=10, color='green')
        axs[0, 0].plot([60, 80], [150, 150], 'b-', lw=3)
    plt.savefig(f'Out_Img/Control_Sample/galaxy_{galaxy}.svg', format='svg', dpi=1200)
    plt.close()
    return

Tabla = []
header_names = ['Group', 'Gal', 'ID', 'ra', 'dec', 'type', 'CHI2NU']
all_filters = ['g', 'r', 'i', 'z']

for filt in all_filters:
    header_names.append(f'XC_{filt}')
    header_names.append(f'e_XC_{filt}')
    header_names.append(f'YC_{filt}')
    header_names.append(f'e_YC_{filt}')
    header_names.append(f'RE_{filt}')
    header_names.append(f'e_RE_{filt}')
    header_names.append(f'MAG_{filt}')
    header_names.append(f'e_MAG_{filt}')
    header_names.append(f'n_{filt}')
    header_names.append(f'e_n_{filt}')
    header_names.append(f'AR_{filt}')
    header_names.append(f'e_AR_{filt}')
    header_names.append(f'PA_{filt}')
    header_names.append(f'e_PA_{filt}')
    
L_galaxies = Table.read('/home/seba/Documents/MorphoLS/Catalog/Galaxies_DECALS_control_sample.csv')
no_ajustados = pd.read_csv('/home/seba/Documents/numeros_unicos.txt', header=None)
n = no_ajustados[0].to_list()
# Filtrar los elementos que NO están en la lista n
L_filt = L_galaxies[~np.isin(L_galaxies['index'], n)]

# Filtrar los elementos que son menores o iguales a 740
L = L_filt[L_filt['index'] <= 740]
print(L)
Datos_L = L.group_by('index')
Grupos = Datos_L.groups.keys
failed_fits=[]
for g in Grupos['index']:
    if g <= 740:
        mask = Datos_L.groups.keys['index'] == g
        Tablef = Datos_L.groups[mask]
        fil_name = filter_sel(g)[0]
        n_filtros = len(fil_name)
        #print(fil_name)
    
        
        # Leer el archivo FITS correspondiente al grupo
        fi = fits.open(f'galfitm_output/Control_Sample/galfitm_galaxy_{g}.fits')
        for i in range(len(Tablef)):
            gal = Tablef[i]['index']
            ID = Tablef[i]['objid']
            ra = Tablef[i]['ra']
            dec = Tablef[i]['dec']
            t = Tablef[i]['type']
            # Leer y añadir datos a la tabla
            header_data = leer_header(i + 1, fi[4].header, fil_name, n_filtros)
            header_data['Group'] = g
            header_data['Gal'] = gal
            header_data['ID'] = ID
            header_data['ra'] = ra
            header_data['dec'] = dec
            header_data['type'] = t
            Tabla.append(header_data)
            # Generar el gráfico si es necesario
            #grafico(fi, g, n_filtros)
    
# Crear la tabla final con los datos y los nombres de los headers
table = Table(rows=Tabla, names=header_names)
ascii.write(table, 'Output_Catalogs/GalfitM_DECALS_Control_Sample_740.csv', format='csv', overwrite=True, fast_writer=False)
print(f'Los grupos no ajustados son {failed_fits}')
print(len(failed_fits))










