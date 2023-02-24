import numpy as np
from numba import njit

@njit
def run_population(tmax=300.0, sampling_time=1.0, cells=int(1000.0), kr=5, kp=10, gamma_r=0.2, gamma_p=0.4, ) -> list:
    cells_arr = []
    for cell_indx in range(1, cells + 1):
        species = np.array([0.0, 0.0, 0.0, 0.0], dtype=np.float64)
        species[1] = cell_indx
        # Constants
        # Reaction matrix
        reaction_type = np.array([[0, 0, 1.0, 0], [0, 0, 0, 1.0], [0, 0, -1.0, 0], [0, 0, 0, -1.0]], dtype=np.int64)
        # Propensities initiation
        propensities = np.zeros(4, dtype=np.float64)
        tarr = np.arange(0, tmax,   sampling_time, dtype=np.float64)
        # Simulation Space
        sim  = np.zeros((len(tarr), len(species)), dtype=np.float64)
        sim[0] = species
        for indx_dt in range(1, len(tarr)):
            species = sim[indx_dt - 1]
            while species[0] < tarr[indx_dt]:
                # Species
            
                r = species[2]
                p = species[3]
                # Propensities
            
                propensities[0] = kr
                propensities[1] = kp * r
                propensities[2] = gamma_r * r
                propensities[3] = gamma_p * p
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