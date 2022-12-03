# ToyMaker

This library is designed to facilitate research of stochasticity effects in genetic circuits. 

## Usage:

ToyMaker uses 6 types of parameters: k_values, propensities, species, reactions and simulation_parameters.

### k_values:

This parameter is for setting the constants that you're going to use in your model. 

### propensities:

This  parameter is used to set the propensities calculation equations. 

### species:

This parameter is to set the species that are involved in the system. 

### reactions:

This parameter is to set what kind of reaction occurs (create, destroy).

### simulation_parameters:

This parameter is used to set maximum simulation time (tmax), lenght per step (sampling_time), and population ammount (cells) to run the simulation.

### Setting up the model: 

Once the model is stablished it's time to set up the model by using: 

    model = tmk.ToyModel(toy_name='My_first_model')
    model.setup(species, reactions, propensities, k_values, simulation_parameters)
    
    
 This generates a Pythonscript with your implemented system, take into account that the generated python file name is the toy_name parameter in ToyModel. 
 Finally, you can use your generated file in a Jupyter Notebook.
 
 ### Example: 
 
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
 
