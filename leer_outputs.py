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
from ejecutable import L

#fil_name = np.array(['g', 'r', 'i', 'z'])

def leer_header(NGAL, HEADER):
    data=[]
    data.append(HEADER['CHI2NU'])
    for i in range(n_filtros):
        xc = HEADER.get(f'{NGAL}_XC_{fil_name[i]}', '0.0 0.0 0.0').split( )
        yc = HEADER.get(f'{NGAL}_YC_{fil_name[i]}', '0.0 0.0 0.0').split( )
        R = HEADER.get(f'{NGAL}_RE_{fil_name[i]}', '0.0 0.0 0.0').split( ) 
        m = HEADER.get(f'{NGAL}_MAG_{fil_name[i]}', '0.0 0.0 0.0').split ()
        N = HEADER.get(f'{NGAL}_n_{fil_name[i]}', '0.0 0.0 0.0').split( )
        ar = HEADER.get(f'{NGAL}_AR_{fil_name[i]}', '0.0 0.0 0.0').split( )
        PA = HEADER.get(f'{NGAL}_PA_{fil_name[i]}', '0.0 0.0 0.0').split( )
        '''xc = HEADER[f'{NGAL}_XC_{fil_name[i]}'].split( )
        yc = HEADER[f'{NGAL}_YC_{fil_name[i]}'].split( )
        R = HEADER[f'{NGAL}_RE_{fil_name[i]}'].split( )
        m = HEADER[f'{NGAL}_MAG_{fil_name[i]}'].split( )
        N = HEADER[f'{NGAL}_n_{fil_name[i]}'].split( )
        ar = HEADER[f'{NGAL}_AR_{fil_name[i]}'].split( )
        PA = HEADER[f'{NGAL}_PA_{fil_name[i]}'].split( )'''
        data.append(xc[0])
        data.append(xc[2])
        data.append(yc[0])
        data.append(yc[2])
        data.append(m[0])
        data.append(m[2])
        data.append(R[0])
        data.append(R[2])
        data.append(N[0])
        data.append(N[2])
        data.append(ar[0])
        data.append(ar[2])
        data.append(PA[0])
        data.append(PA[2])
    return(data)

def grafico(fits, grupo, n_filtros):
    fig, axs = plt.subplots(3, 4, figsize=(20, 8), sharex=True, sharey=True)
    for j in range(n_filtros):
        axs[0, j].imshow(fits[j].data, cmap='gray', vmin = -0.5, vmax = 0.5)
        axs[0, j].set_title(f'{fil_name[j]}')
        axs[1, j].imshow(fits[n_filtros+j].data, cmap='gray', vmin = -0.5, vmax = 0.5)
        axs[2, j].imshow(fits[2*n_filtros+j].data, cmap='gray', vmin = -0.5, vmax = 0.5)
        axs[0, 0].text(60, 140, '6\"', fontsize = 10, color='green')
        axs[0, 0].plot([60, 80], [150, 150], 'b-', lw=3)
    plt.savefig(f'Out_Img/group_{grupo}.svg', format='svg', dpi = 1200)
    plt.close()
    return

Tabla=[]
header_names=['CHI2NU']
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
#L_try = Table.read('/home/seba/Documents/DECALS/Galaxies/Galaxies_group_36.csv')
Datos_L = L.group_by('Group')
Grupos = Datos_L.groups.keys
for g in Grupos['Group']:
    mask = Datos_L.groups.keys['Group'] == g
    Tablef = Datos_L.groups[mask]
    fil_name = filter_sel(g)[0]
    #n_filtros = len(filter_sel(g))
    n_filtros = len(fil_name)
    fi = fits.open(f'galfitm_output/galfitm_group_{g}.fits')
    for i in range(len(Tablef)): 
        Tabla.append(leer_header(i+1, fi[4].header))
        #grafico(fi, g, n_filtros)

'''header_names=['CHI2NU']
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
    header_names.append(f'e_PA_{filt}')'''

table = Table(rows=Tabla, names=header_names)
ascii.write(table, f'GalfitM_output_DECALS_SPLUS_17.csv', format='csv', fast_writer=False)












