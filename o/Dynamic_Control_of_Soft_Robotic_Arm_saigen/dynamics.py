import numpy as np
from math import sin, cos, sqrt

from kinematics import KinematicsOfOneSection


def T2(M):
    return M[0,0] + M[1,1]





class Dynamics(KinematicsOfOneSection):
    """ソフトロボットの動力学"""
    
    K = 1700
    D = 110
    m = 0.13
    
    g = -9.81
    
    
    def __init__(self,):
        
        super().__init__()
        return
    
    
    def update(self, q):
        """いろいろ更新"""
        
        self.dRdq_all = [
            [self.linearized_dRdl1(q, xi) for xi in self.xi_all],
            [self.linearized_dRdl2(q, xi) for xi in self.xi_all],
            [self.linearized_dRdl3(q, xi) for xi in self.xi_all],
        ]
        
        self.dpdq_all = [
            [self.linearized_dpdl1(q, xi) for xi in self.xi_all],
            [self.linearized_dpdl2(q, xi) for xi in self.xi_all],
            [self.linearized_dpdl3(q, xi) for xi in self.xi_all],
        ]
        
        
        
        return
    
    
    def M_omega_jk(j, k):
        pass






if __name__ == "__main__":
    print("hello!")