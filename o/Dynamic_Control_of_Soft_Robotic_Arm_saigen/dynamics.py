import numpy as np
from math import sin, cos, sqrt

from kinematics import KinematicsOfOneSection



class Dynamics(KinematicsOfOneSection):
    """ソフトロボットの動力学"""
    
    K = 1700
    D = 110
    m = 0.13
    
    g = -9.81
    Ixx = 1 # ?
    
    def __init__(self,):
        
        super().__init__()
        
        self.Gv = np.array([[0, 0, self.g]]).T
        
        

        
        return
    
    def update_state(self, q):
        """状態を更新"""
        self.q = q
        
        self.dpdq = [
            [self.linearized_dpdl1(self.q, xi) for xi in self.xi_all],
            [self.linearized_dpdl2(self.q, xi) for xi in self.xi_all],
            [self.linearized_dpdl3(self.q, xi) for xi in self.xi_all],
        ]
        
        self.dRdq = [
            [self.linearized_dRdl1(self.q, xi) for xi in self.xi_all],
            [self.linearized_dRdl2(self.q, xi) for xi in self.xi_all],
            [self.linearized_dRdl3(self.q, xi) for xi in self.xi_all],
        ]
    
    
    def M_omega_jk(self, j, k):
        dRdq = [
            T2(self.dRdq[j+1][i] @ self.dRdq[k+1][i]) for i, _ in enumerate(self.xi_all)
        ]
        
        return self.Ixx * np.sum(dRdq)
    
    
    def M_v_jk(self, j, k):
        dpdq = [
            self.dpdq[j+1][i] @ self.dpdq[k+1][i] for i, _ in enumerate(self.xi_all)
        ]
        
        return self.m * np.sum(dpdq)


    def M_jk(self, j, k):
        return self.M_omega_jk(j, k) + self.M_v_jk(j, k)


    def M(self,):
        """慣性行列"""
        
        M_ = np.zeros((3, 3))
        for j in range(3):
            for k in range(3):
                M_[j, k] = self.M_jk(j, k)
        
        return M_




if __name__ == "__main__":
    hoge = Dynamics()
    hoge.update_state(np.array([[1, 2, 3]]).T)
    M = hoge.M()
    print("hello!")