import numpy as np
from numba import njit

@njit
def run_population(tmax={tmax_value}, sampling_time={sampling_time_value}, cells=int({cells_value}), {k_values_to_write}) -> list:
    cells_arr = []
    for cell_indx in range(1, cells + 1):
        species = np.array({species_values}, dtype=np.float64)
        species[1] = cell_indx
        # Constants
        # Reaction matrix
        reaction_type = np.array({reaction_matrix_values}, dtype=np.int64)
        # Propensities initiation
        propensities = np.zeros({reactions_keys_values}, dtype=np.float64)
        tarr = np.arange(0, tmax,   sampling_time, dtype=np.float64)
        # Simulation Space
        sim  = np.zeros((len(tarr), len(species)), dtype=np.float64)
        sim[0] = species
        for indx_dt in range(1, len(tarr)):
            species = sim[indx_dt - 1]
            while species[0] < tarr[indx_dt]:
                # Species
            {species_to_write}
                # Propensities
            {propensities_to_write}
                τarr = np.zeros(len(propensities), dtype=np.float64)
                # Calculate tau times
                for indx_τ in range(len(propensities)):
                    if propensities[indx_τ] > 0 :
                        τarr[indx_τ] = -(1/propensities[indx_τ]) * np.log(np.random.rand())
                    else:
                        τarr[indx_τ] = np.inf
                τ = np.min(τarr)
                q = np.argmin(τarr)
                species = species + reaction_type[q] if τ != np.inf else species
                species[0] = species[0] + τ
            sim[indx_dt] = species
        cells_arr.append(sim)
    return cells_arr