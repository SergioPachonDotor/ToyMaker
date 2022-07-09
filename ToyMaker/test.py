from stochastic_sims.gillespie import cell
from stochastic_sims.stochastic_division import cell_division
from stochastic_sims.tools.stochastic_division_tools import birth_size
import numpy as np

class Toy:

    def __init__(self, name='Chiquimodelo'):
        self.name = name
        self.result = np.array([])
        self.species = {}
        self.reactions = {}
        self.tmax = 100
        self.sampling_time = 1.


    def simulate_gillespie(self):
        self.result = cell(species=self.species, reactions=self.reactions, tmax=self.tmax, sampling_time=self.sampling_time)
        
    
    def simulate_division(self, doubling_time:float=18., noise_at_division:bool=False):
        self.result = cell_division(self.species, self.reactions, self.tmax, self.sampling_time, doubling_time=doubling_time, noise_at_division=noise_at_division)
    

    def make(
            self, 
            species:dict, 
            reactions:dict, 
            sim_type:str='classic', 
            size_b:float=1.0, 
            stochastic_birth_size:bool=False,
            tmax = 100.,
            sampling_time = 1.
            ):

        self.tmax = tmax
        self.sampling_time = sampling_time

        if sim_type == 'classic':
            self.species['time'] = 0.
            self.species['cell'] = 0.
        
        elif sim_type == 'division':
            self.species['time'] = 0.
            self.species['cell'] = 0.
            self.species['size'] = size_b if stochastic_birth_size == False else birth_size()

        for ks,vs in species.items():
            self.species[ks] = vs
        
        for kr,vr in reactions.items():
            self.reactions[kr] = vr
                       

if __name__ == '__main__':
    def k_r():
        return 2

    def gamma_r(Xr):
        return (1/10) * Xr

    species = {
                'Xr': 0.,
    }

    reactions = {

        k_r: {'create': ['Xr']},
        gamma_r:{'destroy': ['Xr']},
        'division': {'segregate': ['Xr']}
    }

    tmax = 40
    sampling_time = 1
    cells = 1
    model = Toy(name='Test')
    model.make(species=species, reactions=reactions, tmax=100, sampling_time=1)