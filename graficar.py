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

def grafico(fits, grupo, n_filtros, crop_size_arcmin=3):
    pix_scale = 0.262  # arcsec/pixel
    crop_size_pix = int((crop_size_arcmin * 60) / pix_scale)  # De arcmin a píxeles
    half_crop = crop_size_pix // 2

    fig, axs = plt.subplots(3, 4, figsize=(20, 8), sharex=True, sharey=True)

    for j in range(n_filtros):
        # Obtener tamaño original de la imagen
        data = fits[j].data
        y_size, x_size = data.shape

        # Definir los límites del recorte
        x_min, x_max = (x_size // 2 - half_crop, x_size // 2 + half_crop)
        y_min, y_max = (y_size // 2 - half_crop, y_size // 2 + half_crop)

        # Aplicar recorte
        cropped_data_1 = fits[j].data[y_min:y_max, x_min:x_max]
        cropped_data_2 = fits[n_filtros + j].data[y_min:y_max, x_min:x_max]
        cropped_data_3 = fits[2 * n_filtros + j].data[y_min:y_max, x_min:x_max]

        # Graficar imágenes recortadas
        axs[0, j].imshow(cropped_data_1, cmap='gray', vmin=-0.5, vmax=0.5)
        axs[0, j].set_title(f'{fil_name[j]}')

        axs[1, j].imshow(cropped_data_2, cmap='gray', vmin=-0.5, vmax=0.5)
        axs[2, j].imshow(cropped_data_3, cmap='gray', vmin=-0.5, vmax=0.5)

        # Añadir barra de escala en la primera imagen
        axs[0, 0].text(30, 60, '1 arcmin', fontsize=10, color='green')
        axs[0, 0].plot([30, 259], [80, 80], color='green', lw=1)
    plt.tight_layout()
    # Reducir los espacios entre las columnas
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    plt.savefig(f'Out_Img/group_{grupo}_cropped.svg', format='svg', dpi=1200)
    plt.close()

Datos_L = L.group_by('Group')
Grupos = Datos_L.groups.keys
for g in Grupos['Group']:
    if g == 235:
        mask = Datos_L.groups.keys['Group'] == g
        Tablef = Datos_L.groups[mask]
        fil_name = filter_sel(g)[0]
        n_filtros = len(fil_name)
        #print(fil_name)
    
        try:
            # Leer el archivo FITS correspondiente al grupo
            fi = fits.open(f'galfitm_output_2/galfitm_group_{g}.fits')
            for i in range(len(Tablef)):
                # Generar el gráfico si es necesario
                grafico(fi, g, n_filtros)
        except Exception as e:
            print(f'El grupo {g} no fue graficado, {e}')

