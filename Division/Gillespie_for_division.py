import numpy as np
from MakerTools.tools import *
from MakerTools.division_tools import *


def Simulate_Division(species:dict, reactions:dict, tmax:int, sampling_time:float=0.1, cell=1, division=18):
    
    mu = np.log(2)/division

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
        
        sim[i], time, division_time = Gillespie_for_division(sim[i - 1], tarr[i], reaction_type, propensities, reactions_species_index, mu, division_time, division_index)
        division_time -= time
        sim[i][2] *= np.exp(mu * time)
        sim[i][1] = cell

    return sim


def Gillespie_for_division(species, tmax, reaction_type, propensities, reactions_species_index, mu, division_time, division_index):
    time = 0.
    while species[0] < tmax:
        # τarr = calculate_propensities(propensities, species, reactions_species_index)
        τarr = calculate_propensities_for_division(propensities, species, reactions_species_index, mu, species[2], division_index)
        τ,  q = minimal_value(τarr)

        if division_time > τ:
            species += reaction_type[q]
            species[2] *= np.exp(mu * τ)
            division_time -= τ

        elif division_time < τ:
            species[2] /= 2
            species = dilute(species, reaction_type, division_index)
            division_time = time_to_division(mu, species[2])

        time = tmax - species[0]
        species[0] += τ

    return species, time, division_time