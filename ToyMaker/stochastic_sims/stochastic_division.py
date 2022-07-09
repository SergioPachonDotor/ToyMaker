from stochastic_sims.tools.stochastic_division_tools import *
from stochastic_sims.tools.common_tools import *


def cell_division(species:dict, reactions:dict, tmax:int, sampling_time:float=0.1, cell:float=1., doubling_time:float=18., noise_at_division=False):
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

        sim[i], division_time = gillespie_for_division(sim[i - 1], tarr[i], reaction_type, propensities, reactions_species_index, mu, division_time, division_index, noise_at_division)
        sim[i][1] = cell
        sim[i][0] = np.float64(sim[i][0])
    return sim


def gillespie_for_division(species, reference_time, reaction_type, propensities, reactions_species_index, mu, division_time, division_index, noise_at_division):
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