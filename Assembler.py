import numpy as np
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
        propensities_to_write   = ''.join([ f'\n\t\t\t\tpropensities[{k}] = {v}'   for k,v in pindex.items()])
        k_values_to_write       = ''.join([ f'\n\t\t{k} = {v}'                 for k,v in self.k_values.items()])
        species_to_write        = ''.join([ f'\n\t\t\t\t{v} = species[{k}]'        for k,v in dict(enumerate(self.species.keys())).items() if k > 1])
        #
        gillespie_str= f"import numpy as np\
                        \nfrom numba import njit\n\
                        \n@njit\
                        \ndef run_population(tmax={self.simulation_parameters['tmax']}, sampling_time={self.simulation_parameters['sampling_time']}, cells={int(self.simulation_parameters['cells'])}) -> np.array:\
                        \n\tcells_arr = []\
                        \n\tfor cell_indx in range(1, cells + 1):\
                        \n\n\t\tspecies = np.array({list(self.species.values())}, dtype=np.float64)\
                        \n\t\tspecies[1] = cell_indx\
                        \n\n\t\t# Constants\
                        \t\t{k_values_to_write}\
                        \n\n\t\t# Reaction matrix\
                        \n\t\treaction_type = np.array({self.reaction_matrix}, dtype=np.int64)\
                        \n\n\t\t# Propensities initiation\
                        \n\t\tpropensities = np.zeros({len(self.reactions.keys())}, dtype=np.float64)\
                        \n\t\ttarr = np.arange(0, tmax,   sampling_time, dtype=np.float64)\
                        \n\n\t\t# Simulation Space\
                        \n\t\tsim  = np.zeros((len(tarr), len(species)), dtype=np.float64)\
                        \n\t\tsim[0] = species\
                        \n\n\t\tfor indx_dt in range(1, len(tarr)):\
                        \n\t\t\tspecies = sim[indx_dt - 1]\
                        \n\n\t\t\twhile species[0] < tarr[indx_dt]:\
                        \n\t\t\t\t# Species\
                        \t\t\t{species_to_write}\
                        \n\n\t\t\t\t# Propensities\
                        \t\t\t{propensities_to_write}\
                        \n\n\t\t\t\tτarr = np.zeros(len(propensities), dtype=np.float64)\
                        \n\n\t\t\t\t# Calculate tau times\
                        \n\t\t\t\tfor indx_τ in range(len(propensities)):\
                        \n\t\t\t\t\tif propensities[indx_τ] > 0 :\
                        \n\t\t\t\t\t\tτarr[indx_τ] = -(1/propensities[indx_τ]) * np.log(np.random.rand())\
                        \n\t\t\t\t\telse:\
                        \n\t\t\t\t\t\tτarr[indx_τ] = np.inf\
                        \n\n\t\t\t\tτ = np.min(τarr)\
                        \n\t\t\t\tq = np.argmin(τarr)\
                        \n\t\t\t\tspecies = species + reaction_type[q] if -1 not in species + reaction_type[q] else species\
                        \n\t\t\t\tspecies[0] = species[0] + τ\
                        \n\t\t\tsim[indx_dt] = species\
                        \n\t\tcells_arr.append(sim)\
                        \n\treturn cells_arr"
            
        return gillespie_str

    def assemble(self) -> None:
        data = self.create_gillespie_string()

        if os.path.exists('./results/'):
            pass

        else:
            os.mkdir('./results/')
            
        with open(f'./results/{self.toy_name}.py', 'w') as f:
            f.writelines(
                f"{data}")

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

    cells = 1
    tmax = 10
    sampling_time = 1

    model = ToyModel(toy_name='My_first_model')
    model.setup(species, reactions, propensities, k_values, simulation_parameters)