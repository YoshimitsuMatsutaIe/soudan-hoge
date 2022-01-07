import numpy as np


import N_is_1.M0_0
import N_is_1.M0_1
import N_is_1.M0_2
import N_is_1.M1_0
import N_is_1.M1_1
import N_is_1.M1_2
import N_is_1.M2_0
import N_is_1.M2_1
import N_is_1.M2_2


def f(q_large, xi_large, q_dot_large):
    """wow"""
    m_0_0 = N_is_1.M0_0.f(q_large, xi_large, q_dot_large)
    m_0_1 = N_is_1.M0_1.f(q_large, xi_large, q_dot_large)
    m_0_2 = N_is_1.M0_2.f(q_large, xi_large, q_dot_large)
    m_1_0 = N_is_1.M1_0.f(q_large, xi_large, q_dot_large)
    m_1_1 = N_is_1.M1_1.f(q_large, xi_large, q_dot_large)
    m_1_2 = N_is_1.M1_2.f(q_large, xi_large, q_dot_large)
    m_2_0 = N_is_1.M2_0.f(q_large, xi_large, q_dot_large)
    m_2_1 = N_is_1.M2_1.f(q_large, xi_large, q_dot_large)
    m_2_2 = N_is_1.M2_2.f(q_large, xi_large, q_dot_large)
    
    return np.array([
        [m_0_0, m_0_1, m_0_2],
        [m_1_0, m_1_1, m_1_2],
        [m_2_0, m_2_1, m_2_2],
    ])


