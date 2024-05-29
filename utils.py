import os
import numpy as np
def filter_sel(grupo):
    filtros_existentes=[]
    lambda_por_filtro=[]
    filtros_longitudes = {'g': 4770, 'r': 6231, 'i': 7625, 'z': 9134}
    directorio = f'/home/seba/Documents/DECALS/joined_bricks/{grupo}/'
    for filtro, lambda_ in filtros_longitudes.items():
        nombre_archivo = f'{grupo}_image_{filtro}.fits'
        if  nombre_archivo in os.listdir(directorio):
            filtros_existentes.append(filtro)
            lambda_por_filtro.append(lambda_)
    filtros_existentes = np.array(filtros_existentes)
    lambda_por_filtro = np.array(lambda_por_filtro)
    return filtros_existentes, lambda_por_filtro

