import numpy as np

import C
import G
import M


q_large = np.array([[0.01, 0.01, 0.02]]).T
q_dot_large = np.array([[0.00, 0.00, 0.00]]).T
xi_large = np.array([[1, 1, 1]]).T


print(G.f(q_large, xi_large, q_dot_large))