import numpy as np

import N_is_1.G0
import N_is_1.G1
import N_is_1.G2



def f(q_large, xi_large, q_dot_large):
    g0 = N_is_1.G0.f(q_large, xi_large, q_dot_large)
    g1 = N_is_1.G1.f(q_large, xi_large, q_dot_large)
    g2 = N_is_1.G2.f(q_large, xi_large, q_dot_large)
    
    return np.array([
        [g0],
        [g1],
        [g2],
    ])


