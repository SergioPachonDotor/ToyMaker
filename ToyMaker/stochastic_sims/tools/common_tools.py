import numpy as np

# Common Simulation Tools
def minimal_value(array):
    """
    Choose the lowest propensity value and returns the 
    tau value and the q value.
    """
    x = array[0]
    indx = 0

    for i in range(1,len(array)):
        if array[i] < x:
            x = array[i]
            indx = i
    return x,  indx


def change_matrix(reaction_model, species) -> np.array:
    """
    Returns a change matrix for each reaction case. 
             1 for "create" statement.
            -1 for "desrtoy" statement.
            100 for "burst" statement.
        
        it's essential to respect statements and be sure taht you're
        writing them in the correct fashion.

        Example of change matrix for the following system:

        species = {
                    't':      0., 
                    'cell':   0.,
                    'size':   1.0,
                    'z':      0.,
                    'r':      0.,
                    'p':      0.,       
                    }

        reactions = {
                        kz: {'create':    ['z']},
                        kr: {'create':    ['r']},
                        kp: {'create':    ['p']},
                        gamma_z : {'destroy':   ['z']},
                        gamma_r : {'destroy':   ['r']},
                        gamma_p : {'destroy':   ['p']},
                        division: {'segregate' : ['z', 'p', 'r']}
                        }

           [[ 0.  0.  0.  1.  0.  0.]
            [ 0.  0.  0.  0.  1.  0.]
            [ 0.  0.  0.  0.  0.  1.]
            [ 0.  0.  0. -1.  0.  0.]
            [ 0.  0.  0.  0. -1.  0.]
            [ 0.  0.  0.  0.  0. -1.]]

    """
    reactions: np.zeros = np.zeros((len(reaction_model.keys()),len(species)))
    
    reaction_type: list = list(reaction_model.values())

    for row in range(len(reactions)):
        for column in range(len(reactions[row])):

            if  'create' in reaction_type[row].keys() and species[column] in reaction_type[row]['create']:
                reactions[row][column] = 1.

            elif 'destroy' in reaction_type[row].keys() and species[column] in reaction_type[row]['destroy']:
                reactions[row][column] = -1.

            elif 'burst' in reaction_type[row].keys() and species[column] in reaction_type[row]['burst']:
                reactions[row][column] = 100.

            elif 'segregate' in reaction_type[row].keys() and species[column] in reaction_type[row]['segregate']:
                reactions[row][column] = 5.

            elif 'create_mrna' in reaction_type[row].keys() and species[column] in reaction_type[row]['create_mrna']:
                reactions[row][column] = 7.
    
    return reactions


def get_funct_get_pams(fun_dict):

    """ Gets the function parameters and it's indexes"""

    function_arr = np.zeros(len(fun_dict))
    
    fun_dict = dict(enumerate(fun_dict.keys())).items()
    function_arr = [fun for idx, fun in fun_dict if type(fun) != str]
    function_idx = [idx for idx, fun in fun_dict if type(fun) != str]
    function_parameters_arr = [s.__code__.co_varnames[:s.__code__.co_argcount] for s in function_arr]
    return function_arr, function_parameters_arr, function_idx


def get_species_names(model):
    return [k if type(k) == str else k.__name__ for k in model.keys()]


def get_index(names, parameters, system):
    return [[names.index(k) for k in parameters[i]] for i in range(len(system))]


def set_values(model) -> dict:
    """Returns a dictionary of type = {Function Name: Default Value}"""
    return {k if type(k) == str else k.__name__ : v  for k,v in model.items()}
