def k_r():
    return 1

def gamma_r():
    return 1/5

species = {
            't':    0., 
            'cell': 0,
            'Xr': 0.,
}
reactions = {

    k_r: {'create': ['Xr']},
    gamma_r:{'destroy': ['Xr']},

}
