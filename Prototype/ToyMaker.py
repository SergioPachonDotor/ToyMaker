import numpy as np
import warnings
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
# ___________________________________________ 
# 
# Code Wrtitten By Sergio Andrés Pachón Dotor
# 2022
# e-mail: sap9827@gmail.com
# 
# ___________________________________________
#
# Inspired in specific Gillespie from: https://github.com/jdmarmolejo/Propagation-of-Noise-in-Feedback-Gene-Networks/blob/main/X%20bloquea%20a%20Y%20y%20Y%20activa%20a%20X/ABloquea%20BActiva.ipynb
# Multiprocessing solution taken from: https://analyticsindiamag.com/run-python-code-in-parallel-using-multiprocessing/

def get_funct_get_pams(fun_dict):
    function_arr = np.zeros(len(fun_dict))
    
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
                        division: {'segregate' : ['z', 'p', 'r']}
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

            elif 'segregate' in reaction_type[row].keys() and species[column] in reaction_type[row]['segregate']:
                reactions[row][column] = 5

            elif 'create_mrna' in reaction_type[row].keys() and species[column] in reaction_type[row]['create_mrna']:
                reactions[row][column] = 7
    
    return reactions


def setup_sim(tarr, species, cell):
    """
    Returns a time array of the simulation with 
    default values already set.
    """
    sim = np.zeros((len(tarr),len(species)), dtype=np.float64)
    init_state = np.array(list(species.values()), dtype=np.float64)
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
    x = array[0]
    indx = 0
    for i in range(1,len(array)):
        if array[i] < x:
            x = array[i]
            indx = i
    return x,  indx


def calculate_tau(p):
    return -(1/p) * np.log(np.random.rand()) if p > 0. else np.inf


def calculate_tau_for_division(p, mu, size):
    return (1/mu) * np.log(1 - (mu * np.log(np.random.rand())/(p * size))) if p > 0. else np.inf


def calculate_propensities(propensities, species, reactions_species_index):
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


def calculate_propensities_for_division(propensities, species, reactions_species_index, mu, size, division_index):
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


def birth_size(size=1.0) -> float:
    return np.random.normal(loc=size, scale=0.1)


def grow(mu, time, size):
    """
        Makes cell grow at rate s0 * e^(mu * t)
    """
    return size * np.exp(mu * time)


def divide(size, noise_at_division=False) -> tuple:
    """Returns Size, β"""
    if noise_at_division == False:
        β = 0.5
        size = size / 2.
        return size, β

    elif noise_at_division == True:
        β = np.random.beta(50, 50)
        size = size * β
        return  size, β


def segregate(species, reaction_type, division_index, beta):
    """Segregation"""
    for i in range(len(species)):
        if reaction_type[division_index][i] == 5:
            species[i] = np.random.binomial(species[i], beta)
    return species


# Classic Gillespie

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
        sim[i][1] = cell
        
    return sim


def Gillespie(species, tmax, reaction_type, propensities, system, species_index, reactions_species_index, system_idx, deterministic):
    while species[0] < tmax:
        τarr = calculate_propensities(propensities, species, reactions_species_index)
        τ,  q = minimal_value(τarr)
        species = species + reaction_type[q]
        species[0] = species[0] + τ
    return species


# Gillespie for Division Studies

def Simulate_Division(species:dict, reactions:dict, tmax:int, sampling_time:float=0.1, cell=1, doubling_time=18, noise_at_division=False):
    # Size params setup

    mu = np.log(2)/doubling_time

    species_names = get_species_names(species)
    division_index = get_species_names(reactions).index('division')

    reaction_type = change_matrix(reactions, species_names)

    propensities, propensities_parameters, params_idx = get_funct_get_pams(reactions)

    species:dict = set_values(species)
    reactions:dict = set_values(reactions)

    reactions_species_index = get_index(species_names, propensities_parameters, propensities)

    tarr = np.arange(0., tmax , sampling_time ,dtype=np.float64)
    sim = setup_division_sim(tarr, species, cell)

    division_time = time_to_division(mu, sim[0][2])

    for i in range(1,len(tarr)):

        sim[i], division_time = Gillespie_for_division(sim[i - 1], tarr[i], reaction_type, propensities, reactions_species_index, mu, division_time, division_index, noise_at_division)
        sim[i][1] = cell
        sim[i][0] = np.float64(sim[i][0])
    return sim


def Gillespie_for_division(species, reference_time, reaction_type, propensities, reactions_species_index, mu, division_time, division_index, noise_at_division):
    # τarr = calculate_propensities(propensities, species, reactions_species_index)
    time = 0.
    birth_s = species[2]

    τarr = calculate_propensities_for_division(propensities, species, reactions_species_index, mu, species[2], division_index)
    τ,  q = minimal_value(τarr)

    while species[0] + τ < reference_time:     

        if division_time > τ:
            species += reaction_type[q]         
            species[0] = (species[0] + τ) if (species[0] + τ) < reference_time else species[0]
            species[2] = grow(mu, τ, species[2])
            division_time -= τ

        elif division_time < τ:
            # species = segregate(species, reaction_type, division_index)

            species[2], β = divide(species[2], noise_at_division)
            species = segregate(species, reaction_type, division_index, β)
            birth_s = species[2]
            species[2] = grow(mu, division_time, species[2])
            
            species[0] = (species[0] + division_time) if (species[0] + division_time) < reference_time else species[0]
            
            division_time = time_to_division(mu, birth_s)

        τarr = calculate_propensities_for_division(propensities, species, reactions_species_index, mu, species[2], division_index)
        τ,  q = minimal_value(τarr)

    time = reference_time - species[0]
    species[2] = grow(mu, time, species[2])
    division_time = (division_time - time) if (division_time - time) > 0. else 0.
    species[0] = reference_time

    return species, division_time


# Multiple Cells Simulation Types

def multiple_Cells(species, reactions, tmax, sampling_time, cells=10, deterministic=False, file_name='test.csv'):
    
    set_local_to_save(file_name, species)
    for c in tqdm(range(cells)):
        cell = Cell(species, reactions, tmax, sampling_time, cell=c, deterministic=deterministic)
        simple_save(cell, file_name=file_name)


def multiple_sims_for_division(species, reactions, tmax, sampling_time, cells=10, doubling_time=18., file_name='division_test.csv'):
    
    set_local_to_save(file_name, species)
    for c in tqdm(range(cells)):
        cell = Simulate_Division(species, reactions, tmax, sampling_time, cell=c, doubling_time=doubling_time)
        simple_save(cell, file_name=file_name)


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


def get_mean_per_cell_division(cells=[], samples=100, tmax=100, species_idx=3, size_idx=2, single_value=False):

    if single_value == False:
        mean = np.array([np.mean([cells[sample][time:time + 1, species_idx]/cells[sample][time:time + 1, size_idx] for sample in range(samples)]) for time in range(tmax)])
        return mean

    elif single_value == True:
        mean = np.array([np.mean([cells[sample][time:time + 1, species_idx]/cells[sample][time:time + 1, size_idx] for sample in range(samples)]) for time in range(tmax)])
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
    
    def k_r():
        return 2

    def gamma_r(Xr):
        return (1/10)*Xr

    species = {
                't':    0., 
                'cell': 0,
                'size': birth_size(),
                'Xr': 0.,
    }
    reactions = {

        k_r: {'create': ['Xr']},
        gamma_r:{'destroy': ['Xr']},
        'division': {'dilute': ['Xr']}

    }

    tmax = 40
    sampling_time = 1
    cell = 1
    sim = Simulate_Division(species, reactions, tmax, sampling_time, cell, doubling_time=18.)
    # [print(sim[i][0]) for i in range(len(sim))];
    # cell_sim = Cell(species, reactions, tmax, sampling_time, cell=1, deterministic=False)