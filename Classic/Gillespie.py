import numpy as np
from MakerTools.tools import *


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


def gillespie(species, tmax, reaction_type, propensities, reactions_species_index):
    """
    
    """
    while species[0] < tmax:
        τarr = calculate_propensities(propensities, species, reactions_species_index)
        τ,  q = minimal_value(τarr)
        species = species + reaction_type[q]
        species[0] = species[0] + τ
    return species


def simulate_cell(species, reactions, tmax, sampling_time, cell=1):

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
        sim[i] = gillespie(sim[i - 1], tarr[i], reaction_type, propensities, system, species_index, reactions_species_index, system_idx)
        sim[i][0] = round(sim[i][0], 5)
        sim[i][1] = cell
        
    return sim

if __name__ == '__main__':

    def k_r():
        return 1

    def gamma_r():
        return 1/5

    species = {
                't':    0., 
                'cell': 0,
                'Xr': 0.,
    }
    reactions = {

        k_r: {'create': ['Xr']},
        gamma_r:{'destroy': ['Xr']},

    }


    tmax = 100
    sampling_time = 0.1
    cell = 1
    cell_sim = simulate_cell(species, reactions, tmax, sampling_time, cell=1)