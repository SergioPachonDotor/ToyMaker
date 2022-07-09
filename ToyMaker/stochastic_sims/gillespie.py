from stochastic_sims.tools.gillespie_tools import *
from stochastic_sims.tools.common_tools import *

def cell(species, reactions, tmax, sampling_time, cell:int=1):

    species_names = get_species_names(species)

    reaction_type = change_matrix(reactions, species_names)

    propensities, propensities_parameters, params_idx = get_funct_get_pams(reactions)

    species = set_values(species)
    reactions = set_values(reactions)

    reactions_species_index = get_index(species_names, propensities_parameters, propensities)

    tarr = np.arange(0, tmax , sampling_time ,dtype=float) 
    sim = setup_sim(tarr, species, cell)

    for i in range(1,len(tarr)):
        sim[i] = gillespie(sim[i - 1], tarr[i], reaction_type, propensities, reactions_species_index)
        sim[i][1] = cell
        
    return sim


def gillespie(species, tmax, reaction_type, propensities, reactions_species_index):

    while species[0] < tmax:
        τarr = calculate_propensities(propensities, species, reactions_species_index)
        τ,  q = minimal_value(τarr)
        species = species + reaction_type[q]
        species[0] = species[0] + τ
    return species