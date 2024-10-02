from astropy.table import Table
import pandas as pd

cambio = Table.read('/home/seba/Downloads/closest_galaxies_distances_corrected.csv', format='ascii')
L = Table.read('/home/seba/Documents/MorphoLS/Catalogs/Galaxies_Control_Sample.csv', format='ascii')
ajustar = pd.read_csv('/home/seba/Documents/numeros_out.txt', header=None)
n = ajustar[0].to_list()

for galaxia in n:

