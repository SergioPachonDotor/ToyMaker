from ToyMaker import *

# Constants
n1 = 3.5
n2 = 3.5
k1 = 4400 
k2 = 8080
γr1 = 1/5
γr2 = 50/251
β1 = 216/125
β2 = 216/125
kp1 = 50 
kp2 = 45
γp1 = 1/45
γp2 = 1/45

# Species
def p1(p1): return p1
def p2(p2): return p2

# Propensities
def s1 (p2) : return β1 * ((p2**n1) / (k1**n1 + p2**n1))
def s2 (r1) : return γr1 * r1
def s3 (p1) : return β2/(1 + (p1/k2)**n2)
def s4 (r2) : return γr2 * r2
def s5 (r1) : return kp1 * r1
def s6 (p1) : return γp1 * p1
def s7 (r2) : return kp2 * r2
def s8 (p2) : return γp2 * p2

# Species
species = {
            't':    0., 
            'cell': 0, 
            p1:   0., 
            p2:   0., 
            'r1':   0., 
            'r2':   0.}

reactions = {
            s1: {'create':    ['r1']},
            s2: {'destroy':   ['r1']},
            s3: {'create':    ['r2']},
            s4: {'destroy':   ['r2']},
            s5: {'create':    ['p1']},
            s6: {'destroy':   ['p1']},
            s7: {'create':    ['p2']},
            s8: {'destroy':   ['p2']}
}

tmax = 300
sampling_time = 1
samples = 50

# multiple_cells_by_batch(species, reactions, tmax, sampling_time, samples=1, batch_size=10)
multiple_Cells(species, reactions, tmax, sampling_time, cells=samples, file_name='./simulations/test');