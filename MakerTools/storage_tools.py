from tools import *
import os
import pandas as pd

def set_local_to_save(filename, species):
    if os.path.exists(f'{filename}.csv'):
        os.remove(f'{filename}.csv')
        with open(f'{filename}.csv', 'a') as f: 
            header = ','.join(get_species_names(species)) + '\n'
            f.write(header)
    else:
        with open(f'{filename}.csv', 'a') as f: 
            header = ','.join(get_species_names(species)) + '\n'
            f.write(header)


def save_data(cells_arr,file_name):
    
    cells_arr = np.array(cells_arr)
    cells_arr = np.concatenate(cells_arr, axis=0)
    cells_arr = pd.DataFrame(cells_arr)
    
    with open(f'{file_name}.csv', 'a') as f:
        cells_arr.to_csv(f, index=False, header=False, mode='a')


def simple_save(cell,file_name):
    
    cell = np.array(cell)
    cell = pd.DataFrame(cell)
    
    with open(f'{file_name}.csv', 'a') as f:
        cell.to_csv(f, index=False, header=False, mode='a')