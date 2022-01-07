import numpy as np

import C0_0
import C0_1
import C0_2
import C1_0
import C1_1
import C1_2
import C2_0
import C2_1
import C2_2


def f(q_large, xi_large, q_dot_large):
    m_0_0 = C0_0.f(q_large, xi_large, q_dot_large)
    m_0_1 = C0_1.f(q_large, xi_large, q_dot_large)
    m_0_2 = C0_2.f(q_large, xi_large, q_dot_large)
    m_1_0 = C1_0.f(q_large, xi_large, q_dot_large)
    m_1_1 = C1_1.f(q_large, xi_large, q_dot_large)
    m_1_2 = C1_2.f(q_large, xi_large, q_dot_large)
    m_2_0 = C2_0.f(q_large, xi_large, q_dot_large)
    m_2_1 = C2_1.f(q_large, xi_large, q_dot_large)
    m_2_2 = C2_2.f(q_large, xi_large, q_dot_large)
    
    return np.array([
        [m_0_0, m_0_1, m_0_2],
        [m_1_0, m_1_1, m_1_2],
        [m_2_0, m_2_1, m_2_2],
    ])


