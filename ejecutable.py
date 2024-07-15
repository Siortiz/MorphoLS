from astropy.table import *
from astropy.io import ascii
import numpy as np
import csv

L = Table.read('/home/seba/Documents/MorphoLS/Catalog/GalfitM_5_corr.csv')
Datos_L = L.group_by('Group')
Groups = Datos_L.groups.keys

Data = []
#Data = ['python detection.py']
#Data.append('chmod 777 sex_seg.sh')
#Data.append('./sex_seg.sh')
#Data.append('python cross_sex_decals.py')
#Data.append('python psf_mask.py')
#Data.append('python inputs_galfitm.py')

def ejecutable(Groups):
    for g in Groups['Group']:
        Data.append(f'chmod 777 inputs/galfit_{g}.input')
    #    Data.append(f'chmod 777 galfitm-1.4.4-linux-x86_64')
        Data.append(f'./galfitm-1.4.4-linux-x86_64 inputs/galfit_{g}.input')
    #Data.append('python leer_outputs.py')
    fic = open('ejecutable.sh', 'w')
    for line in Data:
        print(line, file=fic)
    fic.close()
    return

ejecutable(Groups)
