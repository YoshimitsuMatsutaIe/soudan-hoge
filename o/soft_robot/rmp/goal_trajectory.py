"""ここに目標軌道のクラスを書く

軌道追従で使用
"""

import numpy as np
from math import pi, cos, sin


def rotate_3d(alpha, beta, gamma):
    """3次元回転行列"""
    
    Rx = np.array([
        [1, 0, 0],
        [0, cos(alpha), sin(alpha)],
        [0, sin(alpha), np.cos(alpha)],
    ])
    Ry = np.array([
        [cos(beta), 0, sin(beta)],
        [0, 1, 0],
        [sin(beta), 0, cos(beta)],
    ])
    Rz = np.array([
        [cos(gamma), -sin(gamma), 0],
        [sin(gamma), cos(gamma), 0],
        [0, 0, 1],
    ])
    
    return Rx @ Ry @ Rz



class Point:
    """点"""
    
    name = 'point'
    
    def __init__(self, **kwargs):
        self.center = np.array(kwargs.pop('center')).reshape(len(kwargs.pop('center')), 1)

    def xd(self, t):
        return self.center
    
    def xd_dot(self, t):
        return np.zeros_like(self.center)

    def xd_dot_dot(self, t):
        return np.zeros_like(self.center)



class Circle:
    """円軌道"""
    
    name = 'circle'
    
    def __init__(self, **kwargs):
        self.r = kwargs.pop('r')
        self.omega = kwargs.pop('omega')
        self.center = np.array(kwargs.pop('center')).reshape(kwargs.pop('center'), 1)
        
        alpha = kwargs.pop('alpha')
        beta = kwargs.pop('beta')
        gamma = kwargs.pop('gamma')
        
        self.R = rotate_3d(alpha, beta, gamma)
    
    
    def xd(self, t):
        return self.center +\
            self.R @ \
                np.array([
                    [r * cos(self.omega * t)],
                    [r * sin(self.omega * t)],
                    [0],
                ])


    def xd_dot(self, t):
        return self.center +\
            self.R @ \
                np.array([
                    [r * -self.omega * sin(self.omega * t)],
                    [r * -self.omega * cos(self.omega * t)],
                    [0],
                ])


    def xd_dot_dot(self, t):
        return self.center +\
            self.R @ \
                np.array([
                    [r * -self.omega**2 * cos(self.omega * t)],
                    [r * -self.omega**2 * sin(self.omega * t)],
                    [0],
                ])


class RoseCurve:
    """バラ曲線"""
    
    name = 'rose_curve'
    
    def __init__(self, **kwargs):
        self.r = kwargs.pop('r')
        self.zeta = kwargs.pop('zeta')
        self.omega = kwargs.pop('omega')
        self.center = np.array(kwargs.pop('center')).reshape(kwargs.pop('center'), 1)

        alpha = kwargs.pop('alpha')
        beta = kwargs.pop('beta')
        gamma = kwargs.pop('gamma')
        
        self.R = rotate_3d(alpha, beta, gamma)

    
    def xd(self, t):
        return self.center + \
            self.R @ \
                np.array([
                    [r * cos(self.zeta*t) * cos(self.omega*t)],
                    [r * cos(self.zeta*t) * sin(self.omega*t)],
                    [0],
                ])


    def xd_dot(self, t):
        return self.center + \
            rotate_3d(self.alpha, self.beta, self.gamma) @ \
                np.array([
                    [r * cos(self.zeta*t) * cos(self.omega*t)],
                    [r * cos(self.zeta*t) * sin(self.omega*t)],
                    [0],
                ])


    def xd(self, t):
        return self.center + \
            rotate_3d(self.alpha, self.beta, self.gamma) @ \
                np.array([
                    [r * cos(self.zeta*t) * cos(self.omega*t)],
                    [r * cos(self.zeta*t) * sin(self.omega*t)],
                    [0],
                ])
