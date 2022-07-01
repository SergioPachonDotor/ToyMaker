import numpy as np


def get_mean_per_cell(cells=[], samples=100, tmax=100, species_idx=2, single_value=False):
    
    if single_value == False:
        return np.array([np.mean([cells[sample][time:time + 1, species_idx] for sample in range(samples)]) for time in range(tmax)])
    
    elif single_value == True:
        mean = np.array([np.mean([cells[sample][time:time + 1, species_idx] for sample in range(samples)]) for time in range(tmax)])
        return mean[-100:].mean()


def get_variance_per_cell(cells=[], samples=100, tmax=100, species_idx=2, single_value=False):
    
    if single_value == False:
        return np.array([np.var([cells[sample][time:time + 1, species_idx] for sample in range(samples)]) for time in range(tmax)])

    elif single_value == True:
        var = np.array([np.var([cells[sample][time:time + 1, species_idx] for sample in range(samples)]) for time in range(tmax)])
        return var[-100:].mean()


def calculate_noise(mean, var):
    return var / pow(mean, 2)