from astropy.table import *
from astropy.io import ascii
import numpy as np
import csv
import pandas as pd

L = Table.read('/home/seba/Documents/MorphoLS/Catalog/Galaxies_DECALS_control_sample_low.csv')
Datos_L = L.group_by('index')
Galaxies = Datos_L.groups.keys
#no_ajustar = pd.read_csv('/home/seba/Documents/numeros_unicos.txt', header=None)
Data = []
#Data = ['python detection.py']
#Data.append('chmod 777 sex_seg.sh')
#Data.append('./sex_seg.sh')
#Data.append('python cross_sex_decals.py')
#Data.append('python psf_mask.py')
#Data.append('python inputs_galfitm_control.py')
#n= no_ajustar[0].to_list()
def ejecutable(Galaxies):
    for g in Galaxies['index']:
        Data.append(f'chmod 777 inputs/Control_Sample/galfit_{g}.input')
    #    Data.append(f'chmod 777 galfitm-1.4.4-linux-x86_64')
        Data.append(f'./galfitm-1.4.4-linux-x86_64 inputs/Control_Sample/galfit_{g}.input')
    #Data.append('python leer_outputs.py')
    fic = open('ejecutable.sh', 'w')
    for line in Data:
        print(line, file=fic)
    fic.close()
    return

ejecutable(Galaxies)
