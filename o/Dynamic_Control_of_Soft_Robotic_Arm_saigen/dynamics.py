import numpy as np
from math import sin, cos, sqrt

from kinematics import Kinematics


class Dynamics(Kinematics):
    """ソフトロボットの動力学"""
    
    K = 1700
    D = 110
    m = 0.13
    
    g = -9.81
    
    
    def __init__(self,):
        
        super().__init__()
        return
    
    
    def update_state(self, q):
        """アクチュエータベクトルを更新"""
        self.q = q
        return
    
    
    def M_omega(j, k):
        
