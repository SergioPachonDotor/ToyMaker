from ast import Break
import numpy as np
import warnings
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import copy

# Author: Sergio Andrés Pachón Dotor
# Inspired in specific Gillespie from: https://github.com/jdmarmolejo/Propagation-of-Noise-in-Feedback-Gene-Networks/blob/main/X%20bloquea%20a%20Y%20y%20Y%20activa%20a%20X/ABloquea%20BActiva.ipynb
# Multiprocessing solution taken from: https://analyticsindiamag.com/run-python-code-in-parallel-using-multiprocessing/

def get_funct_get_pams(fun_dict):
    fun_dict = dict(enumerate(fun_dict.keys())).items()
    function_arr = [fun for idx, fun in fun_dict if type(fun) != str]
    function_idx = [idx for idx, fun in fun_dict if type(fun) != str]
    function_parameters_arr = [s.__code__.co_varnames[:s.__code__.co_argcount] for s in function_arr]
    return function_arr, function_parameters_arr, function_idx


def get_species_names(model):
    return [k if type(k) == str else k.__name__ for k in model.keys()]


def get_index(names, parameters, system):
    return [[names.index(k) for k in parameters[i]] for i in range(len(system))]


def set_values(model) -> dict:
    """Returns a dictionary of type = {Function Name: Default Value}"""
    return {k if type(k) == str else k.__name__ : v  for k,v in model.items()}


def change_matrix(reaction_model, species) -> np.array:
    """
    Returns a change matrix for each reaction case. 
             1 for "create" statement.
            -1 for "desrtoy" statement.
            100 for "burst" statement.
        
        it's essential to respect statements and be sure taht you're
        writing them in the correct fashion.

        Example of change matrix for the following system:

        species = {
                    't':      0., 
                    'cell':   0.,
                    'size':   1.0,
                    'z':      0.,
                    'r':      0.,
                    'p':      0.,       
                    }

        reactions = {
                        kz: {'create':    ['z']},
                        kr: {'create':    ['r']},
                        kp: {'create':    ['p']},
                        gamma_z : {'destroy':   ['z']},
                        gamma_r : {'destroy':   ['r']},
                        gamma_p : {'destroy':   ['p']},
                        division: {'divide' : ['z', 'p', 'r']}
                        }

           [[ 0.  0.  0.  1.  0.  0.]
            [ 0.  0.  0.  0.  1.  0.]
            [ 0.  0.  0.  0.  0.  1.]
            [ 0.  0.  0. -1.  0.  0.]
            [ 0.  0.  0.  0. -1.  0.]
            [ 0.  0.  0.  0.  0. -1.]]

    """
    reactions: np.zeros = np.zeros((len(reaction_model.keys()),len(species)))
    
    reaction_type: list = list(reaction_model.values())

    for row in range(len(reactions)):
        for column in range(len(reactions[row])):

            if  'create' in reaction_type[row].keys() and species[column] in reaction_type[row]['create']:
                reactions[row][column] = 1

            elif 'destroy' in reaction_type[row].keys() and species[column] in reaction_type[row]['destroy']:
                reactions[row][column] = -1

            elif 'burst' in reaction_type[row].keys() and species[column] in reaction_type[row]['burst']:
                reactions[row][column] = 100.

            elif 'dilute' in reaction_type[row].keys() and species[column] in reaction_type[row]['dilute']:
                reactions[row][column] = 2.
    
    return reactions


def setup_sim(tarr, species, cell):
    """
    Returns a time array of the simulation with 
    default values already set.
    """
    sim = np.zeros((len(tarr),len(species)))
    init_state = np.array(list(species.values()))
    init_state[1] = cell
    sim[0] = init_state
    return sim


def setup_division_sim(tarr, species, cell) -> np.array:
    """
    Returns a time array of the simulation with 
    default values already set.
    """
    sim = np.zeros((len(tarr),len(species)))
    init_state = np.array(list(species.values()))
    init_state[1] = cell
    sim[0] = init_state

    return sim


def minimal_value(array):
    """
    Choose the lowest propensity value and returns the 
    tau value and the q value.
    """
    τ = array[0]
    q = 0
    for i in range(1,len(array)):
        if array[i] < τ:
            τ = array[i]
            q = i
    return τ,  q


def calculate_tau(p):
    return -(1/p) * np.log(np.random.rand()) if p > 0 else np.inf


def calculate_tau_for_division(p, mu, size):
    return (1/mu) * np.log(1 - (mu * np.log(np.random.rand())/(p * size)))


def calculate_propensities(propensities, species, reactions_species_index, mode='classic'):
    """ 
    Calculates the propensity of each reaction
    Args:
        propensities (list): List of functions
        species (list): List of species
        species_index (list): List of species indexes
    Returns:
        tauarr (list): List of propensities
    """
    species_values = [species[reactions_species_index[i]] for i in range(len(reactions_species_index))]
    pt = [propensities[i](*species_values[i]) for i in range(len(propensities))]
    τarr = [calculate_tau(p) if p > 0 else np.inf for p in pt]
    return τarr


def calculate_propensities_for_division(propensities, species, reactions_species_index, mu, size):
        species_values = [species[reactions_species_index[i]] for i in range(len(reactions_species_index))]
        pt = [propensities[i](*species_values[i]) for i in range(len(propensities))]
        τarr = [calculate_tau_for_division(p, mu, size) if p > 0 else np.inf for p in pt]
        return τarr 


# system, propensities, species, species_index
def solve_deterministic_model(system, species, species_index, system_idx):
    
    """ 
        Model is a list of functions that are going to be solved
        Species is a list of species names 
    """

    species_values = [species[species_index[i]] for i in range(len(species_index))]
    solved = [system[i](*species_values[i]) for i in range(len(system))]
    for i in range(len(system_idx)):
        species[system_idx[i]] = solved[i]

    return solved


def time_to_division(mu, birth_size):
    return (1/mu) * np.log(((birth_size + np.random.normal(loc=1, scale=0.05))/birth_size))


def update_size(mu:float, time:float, tmax:float) -> float:
    return np.exp(mu * (tmax - time))


def birth_size() -> float:
    return np.random.normal(loc=1.0, scale=0.1)


def dilute(species, reaction_type, division_index):
    return [float(np.random.binomial(species[i], 0.5)) if reaction_type[division_index][i] == 2. else species[i] for i in range(len(species))]


def Gillespie(species, tmax, reaction_type, propensities, system, species_index, reactions_species_index, system_idx, deterministic):
    while species[0] < tmax:
        τarr = calculate_propensities(propensities, species, reactions_species_index)
        τ,  q = minimal_value(τarr)
        species = species + reaction_type[q]
        species[0] = species[0] + τ
    return species


def Gillespie_for_division(species, tmax, reaction_type, propensities, reactions_species_index, mu, division_time, division_index, life_time):
    # Get index of dilution reaction for species dilution
    # reference_time = tmax
    τ = 0.0
    while species[0] < tmax:
        
        τarr = calculate_propensities_for_division(propensities, species, reactions_species_index, mu, species[2])
        τ,  q = minimal_value(τarr)

        if division_time > τ:
            species += reaction_type[q]
            species[2] *= np.exp(mu * τ)
            division_time -= τ
            species[0] += τ

        elif division_time < τ:
            species[2] /= 2
            species = dilute(species, reaction_type, division_index)
            division_time = time_to_division(mu, species[2])
            species[0] += τ

    return species, τ, division_time


def Cell(species, reactions, tmax, sampling_time, cell=1, deterministic=False):

    species_names = get_species_names(species)

    reaction_type = change_matrix(reactions, species_names)

    system, system_parameters, system_idx = get_funct_get_pams(species)
    propensities, propensities_parameters, params_idx = get_funct_get_pams(reactions)

    species = set_values(species)
    reactions = set_values(reactions)

    species_index = get_index(species_names, system_parameters, system)
    reactions_species_index = get_index(species_names, propensities_parameters, propensities)

    tarr = np.arange(0, tmax , sampling_time ,dtype=float) 
    sim = setup_sim(tarr, species, cell)

    for i in range(1,len(tarr)):
        sim[i] = Gillespie(sim[i - 1], tarr[i], reaction_type, propensities, system, species_index, reactions_species_index, system_idx, deterministic)
        sim[i][0] = round(sim[i][0], 5)
        sim[i][1] = cell
        
    return sim
     

def Simulate_Division(species:dict, reactions:dict, tmax:int, sampling_time:float=0.1, cell=1, division_time=18):
    # Size params setup

    mu = np.log(2)/division_time
    life_time = 0.0

    species_names = get_species_names(species)
    division_index = get_species_names(reactions).index('division')

    reaction_type = change_matrix(reactions, species_names)

    propensities, propensities_parameters, params_idx = get_funct_get_pams(reactions)

    species:dict = set_values(species)
    reactions:dict = set_values(reactions)

    reactions_species_index = get_index(species_names, propensities_parameters, propensities)

    tarr = np.arange(0, tmax , sampling_time ,dtype=float) 
    sim = setup_division_sim(tarr, species, cell)

    division_time = time_to_division(mu, sim[0][2])

    for i in range(1,len(tarr)):
        
        sim[i], τ, division_time = Gillespie_for_division(sim[i - 1], tarr[i], reaction_type, propensities, reactions_species_index, mu, division_time, division_index, life_time)
        division_time -= sampling_time
        sim[i][2] *= np.exp(mu * sampling_time)

        sim[i][0] = round(sim[i][0], 5)
        sim[i][1] = cell

    return sim


def multiple_Cells(species, reactions, tmax, sampling_time, cells=10, deterministic=False, file_name='test.csv'):
    
    set_local_to_save(file_name, species)
    for c in tqdm(range(cells)):
        cell = []
        cell = Cell(species, reactions, tmax, sampling_time, cell=c, deterministic=deterministic)
        simple_save(cell, file_name=file_name)
    # cells = [Cell(species, reactions, tmax, sampling_time, cell=c, deterministic=deterministic) for c in tqdm(range(cells))]s
    

def multiple_cells_by_batch(species, reactions, tmax, sampling_time, samples=10, deterministic=False, batch_size=10, file_name='Toy_model'):

    set_local_to_save(file_name, species)
    cells = [i * batch_size for i in range(samples)]

    for i in tqdm(range(samples)):
        cells_arr = []
        with ProcessPoolExecutor() as executor:
            warnings.simplefilter("ignore")
            cells_arr = [executor.submit(Cell, species, reactions, tmax, sampling_time, c + 1, deterministic) for c in range(cells[i], cells[i + 1] if i + 1 < samples else cells[i] + batch_size)]
            cells_arr = [c.result() for c in cells_arr]
            cells_arr = np.array(cells_arr)
            cells_arr = np.concatenate(cells_arr, axis=0)
            cells_arr = pd.DataFrame(cells_arr)
            with open(f'{file_name}.csv', 'a') as f:
                cells_arr.to_csv(f, index=False, header=False, mode='a')


# Statistical Tools

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


# Storage Tools

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
        
        
if __name__ == '__main__':

    # Propensities
    def gamma_r(): return 1/5
    def gamma_p() : return 1/20
    def gamma_z() : return 1/20

    def kz() : return 3.2 * 1/20
    def kr() : return 1
    def kp() : return 50

    # Species

    def z(z): 
        return kz - gamma_z * z 

    def r(z, r):
        return kr * z * (gamma_z/kz) - gamma_r * r

    def p(r, p):
        return kp * r - gamma_p * p

    # Model Deffinition

    species = {
                't':    0., 
                'cell': 0,
                'size': birth_size(),
                z: 0.,
                r: 0.,
                p: 0.,            
                }

    reactions = {
                kz: {'create':    ['z']},
                kr: {'create':    ['r']},
                kp: {'create':    ['p']},
                'division': {'dilute' : ['z', 'p', 'r']},
                gamma_z : {'destroy':   ['z']},
                gamma_r : {'destroy':   ['r']},
                gamma_p : {'destroy':   ['p']},
                
    }

    import matplotlib.pyplot as plt
    tmax = 500
    sampling_time = 0.1
    samples = 10
    cell = 1
    mu = np.log(2)/18
    # cells = multiple_Cells(species, reactions, tmax, sampling_time, cells=10, deterministic=False, file_name='./simulations/Tactic_model')
    sim = Simulate_Division(species, reactions, tmax, sampling_time, cell)
    plt.plot(sim[::][0], sim[::][3])
    plt.savefig('test.jpg')

    # multiple_Cells(species, reactions, tmax, sampling_time, samples=samples, batch_size=100, file_name='./simulations/Tactic_model')