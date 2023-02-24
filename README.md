# ToyMaker 🤖

This library is designed to facilitate research of stochasticity effects in genetic circuits.

## Installation

ToyMaker can be installed using pip:

        pip install ToyMaker

## Usage:

ToyMaker uses 6 types of parameters: k_values, propensities, species, reactions and simulation_parameters.

- ### k_values:

    This parameter is for setting the constants that you're going to use in your model. 

- ### propensities:

    This  parameter is used to set the propensities calculation equations. 

- ### species:

    This parameter is to set the species that are involved in the system. 

- ### reactions:

    This parameter is to set what kind of reaction occurs (create, destroy).

- ### simulation_parameters:

    This parameter is used to set maximum simulation time (tmax), lenght per step (sampling_time), and population ammount (cells) to run the simulation.

- ### Setting up the model: 

    Once the model is stablished it's time to set up the model by using: 

        model = tmk.ToyModel(toy_name='My_first_model')
        model.setup(species, reactions, propensities, k_values, simulation_parameters)
    
    
    This generates a Pythonscript with your implemented system, take into account that the generated python file name is the toy_name parameter in ToyModel. 
    Finally, you can use your generated file in a Jupyter Notebook.
 
 - ### Example: 

        from ToyMaker.Assembler import ToyModel

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
        }

- ### Running the model

    Once My_first_model is setted up, it's time to create a jupyter notebook to import the dynamic library created. Then you can use the "runpopulation" method

        import My_first_model as model
        import numpy as np
        import matplotlib.pyplot as plt

        simulation = np.array(model.run_population(
            tmax = 300, 
            cells=1000, 
            kr=5, 
            gamma_r=1/5, 
            kp=10, 
            gamma_p=0.4))

    Notice that it's possible to change constants values like kr, gamma_r, kp or gamma_p.

- ### Data analysis
    For data analysis it's suggested to convert output of model.runpopulation() to an np.array() to visualize data by using:

        fig, ax = plt.subplots(figsize = (7, 4))
        plt.plot(simulation[:].T[2,:], lw = 0.1, color='red', alpha=0.2) # Multiple Cells
        plt.plot(np.mean(simulation[:].T[2,:], axis=1), color='black') # Mean
        plt.xlabel('Time (a.u.)')
        plt.ylabel('RNA')

![](https://github.com/SergioPachonDotor/ToyMaker/blob/main/ToyMaker/images/example.png)
