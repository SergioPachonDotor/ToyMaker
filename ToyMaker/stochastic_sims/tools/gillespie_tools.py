import numpy as np


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

    
def calculate_tau(p):
    """Calculates next reaction time (τ)"""
    
    return -(1/p) * np.log(np.random.rand()) if p > 0. else np.inf


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


