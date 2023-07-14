import os

class ToyModel:
    
    def __init__(self, toy_name='ToyModel'):
        self.species:dict               = {'time': 0., 'cell': 0.}
        self.reactions:dict             = {}
        self.propensities:dict          = {}
        self.reaction_matrix:list       = []
        self.k_values:dict              = {}
        self.simulation_parameters:dict = {}
        self.toy_name:str               = toy_name

    def add_change_matrix(self, burst_vaue=25.):
        species = [k for k in self.species.keys()]
        reactions = [[0 for c in range(len(species))] for r in range(len(self.reactions))]
        reaction_type:list[dict] = list(self.reactions.values())
        
        for row in range(len(reactions)):
            for column in range(len(reactions[row])):
                if  'create' in reaction_type[row].keys() and species[column] in reaction_type[row]['create']:
                    reactions[row][column] = 1.

                elif 'destroy' in reaction_type[row].keys() and species[column] in reaction_type[row]['destroy']:
                    reactions[row][column] = -1.

                elif 'burst' in reaction_type[row].keys() and species[column] in reaction_type[row]['burst']:
                    reactions[row][column] = burst_vaue
        
        self.reaction_matrix = list(reactions)
        
    def return_model(self) -> dict:
        return self.__dict__
    
    def setup(self, species, reactions, propensities, k_values, simulation_parameters) -> None:
        
        self.species.update(species)
        self.reactions.update(reactions)
        self.propensities.update(propensities)
        self.k_values.update(k_values)
        self.simulation_parameters.update(simulation_parameters)
        self.add_change_matrix()   
        self.create_gillespie_string()
        self.assemble()   
        
    def create_gillespie_string(self) -> str:
        #
        indx = [i for i in range(len(self.propensities.keys()))]    # Propensities Keys by Index
        vals = [s for s in self.propensities.values()]              # Propensities Values
        pindex = dict(zip(indx, vals))                              # Zip index, values
        #
        propensities_to_write   = ''.join([ f'\n                propensities[{k}] = {v}'   for k,v in pindex.items()])
        k_values_to_write       = ''.join([ f'{k}={v}, '                 for k,v in self.k_values.items()])
        species_to_write        = ''.join([ f'\n                {v} = species[{k}]'        for k,v in dict(enumerate(self.species.keys())).items() if k > 1])
        #
        replacements = {
                            "tmax_value" : self.simulation_parameters['tmax'],
                            "sampling_time_value" : self.simulation_parameters['sampling_time'],
                            "cells_value" : self.simulation_parameters['cells'],
                            "species_values" : list(self.species.values()),
                            "reaction_matrix_values" : self.reaction_matrix,
                            "reactions_keys_values" : len(list(self.reactions.keys())),
                            "k_values_to_write": k_values_to_write,
                            "species_to_write": species_to_write,
                            "propensities_to_write": propensities_to_write,
                        }
        
        gillespie_str = ['import numpy as np\n', 
                         'from numba import njit\n', 
                         '\n', 
                         '@njit\n', 
                         'def run_population(tmax={tmax_value}, sampling_time={sampling_time_value}, cells=int({cells_value}), {k_values_to_write}) -> list:\n', 
                         '    cells_arr = []\n', '    for cell_indx in range(1, cells + 1):\n', 
                         '        species = np.array({species_values}, dtype=np.float64)\n', 
                         '        species[1] = cell_indx\n', '        # Constants\n', 
                         '        # Reaction matrix\n', 
                         '        reaction_type = np.array({reaction_matrix_values}, dtype=np.int64)\n', 
                         '        # Propensities initiation\n', 
                         '        propensities = np.zeros({reactions_keys_values}, dtype=np.float64)\n', 
                         '        tarr = np.arange(0, tmax,   sampling_time, dtype=np.float64)\n', 
                         '        # Simulation Space\n', '        sim  = np.zeros((len(tarr), len(species)), dtype=np.float64)\n', 
                         '        sim[0] = species\n', '        for indx_dt in range(1, len(tarr)):\n', 
                         '            species = sim[indx_dt - 1]\n', '            while species[0] < tarr[indx_dt]:\n', 
                         '                # Species\n', '            {species_to_write}\n', '                # Propensities\n', 
                         '            {propensities_to_write}\n', 
                         '                tauarr = np.zeros(len(propensities), dtype=np.float64)\n', 
                         '                # Calculate tau times\n', 
                         '                for indx_tau in range(len(propensities)):\n', 
                         '                    if propensities[indx_tau] > 0 :\n', 
                         '                        tauarr[indx_tau] = -(1/propensities[indx_tau]) * np.log(np.random.rand())\n', 
                         '                    else:\n', 
                         '                        tauarr[indx_tau] = np.inf\n', 
                         '                tau = np.min(tauarr)\n', 
                         '                q = np.argmin(tauarr)\n', 
                         '                species = species + reaction_type[q] if tau != np.inf else species\n', 
                         '                species[0] = species[0] + tau\n', 
                         '            sim[indx_dt] = species\n', 
                         '        cells_arr.append(sim)\n', 
                         '    return cells_arr']
        
        gillespie_str = ''.join(gillespie_str).format(**replacements)   
        return gillespie_str

    def assemble(self) -> None:
        data = self.create_gillespie_string()
        
        with open(f'./{self.toy_name}.py', 'w') as f:
            f.writelines(f"{data}")
        
if __name__ =='__main__':
    
    k_values = {
                'kr' : 5,
                'kp' : 10,
                'gamma_r' : 1/5,
                'gamma_p' : 0.4
            }

    propensities = {
                    'R1': 'kr',
                    'R2': 'kp * r',
                    'R3': 'gamma_r * r',
                    'R4': 'gamma_p * p',
                }

    species = {
            'r':    0.,
            'p':    0.
    }

    reactions = {

        'R1':  {'create':  ['r']},
        'R2':  {'create':  ['p']},
        'R3':  {'destroy': ['r']},
        'R4':  {'destroy': ['p']}
    }

    simulation_parameters = {
        'tmax': 300.,
        'sampling_time': 1.,
        'cells': 1000.
    }

    model = ToyModel(toy_name='My_first_model')
    model.setup(species, reactions, propensities, k_values, simulation_parameters)