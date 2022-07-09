import numpy as np

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


def calculate_tau_for_division(p, mu, size):
    return (1/mu) * np.log(1 - (mu * np.log(np.random.rand())/(p * size))) if p > 0. else np.inf


def calculate_propensities_for_division(propensities, species, reactions_species_index, mu, size, division_index):
    species_values = [species[reactions_species_index[i]] for i in range(len(reactions_species_index))]
    pt = [propensities[i](*species_values[i]) for i in range(len(propensities))]
    τarr = [calculate_tau_for_division(p, mu, size) if p > 0 else np.inf for p in pt]
    return τarr 


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

