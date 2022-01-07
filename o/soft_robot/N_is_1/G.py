import numpy as np

import G0
import G1
import G2



def f(q_large, xi_large, q_dot_large):
    g0 = G0.f(q_large, xi_large, q_dot_large)
    g1 = G1.f(q_large, xi_large, q_dot_large)
    g2 = G2.f(q_large, xi_large, q_dot_large)
    
    return np.array([
        [g0],
        [g1],
        [g2],
    ])


