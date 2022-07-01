import numpy as np

def time_to_division(mu, birth_size):
    return (1/mu) * np.log(((birth_size + np.random.normal(loc=1, scale=0.05))/birth_size))


def calculate_propensities_for_division(propensities, species, reactions_species_index, mu, size, division_index):
    # print(len(propensities), len(species), len(reactions_species_index))
    
    species_values = [species[reactions_species_index[i]] for i in range(len(reactions_species_index)) if i != division_index or i != []]
    pt = [propensities[i](*species_values[i]) for i in range(len(propensities))]
    τarr = [calculate_tau_for_division(p, mu, size) if p > 0 else np.inf for p in pt]
    
    return τarr 


def calculate_tau_for_division(p, mu, size):
    return (1/mu) * np.log(1 - (mu * np.log(np.random.rand())/(p * size)))


def update_size(mu:float, time:float, tmax:float) -> float:
    return np.exp(mu * (tmax - time))


def birth_size() -> float:
    return np.random.normal(loc=1.0, scale=0.1)


def dilute(species, reaction_type, division_index):
    return [float(np.random.binomial(species[i], 0.5)) if reaction_type[division_index][i] == 2. else species[i] for i in range(len(species))]