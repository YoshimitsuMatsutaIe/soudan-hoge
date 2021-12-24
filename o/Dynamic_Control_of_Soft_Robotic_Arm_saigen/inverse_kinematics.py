import numpy as np
from math import sqrt

from scipy.optimize import fsolve

from kinematics import KinematicsOfOneSection


class InverseKinematics(KinematicsOfOneSection):
    """逆運動学を計算"""

    def __init__(self,):
        super().__init__()
        return
    
    
    def calc_qd_from_pd_using_fsolve(self, pd, xi=1):
        """逆運動学を解く（fsolveで）"""
        
        def f(q, xi, pd):
            return self.linearized_mapping_from_actuator_to_task_p(q, xi) - pd
    
        return fsolve(f, np.zeros((3, 1)),args=(xi, pd))
    
    
    
    
    
    def calc_qd_from_pd_1(self, xd, xi=1,):
        """逆運動学を解く"""
        
        
        trial = 1000
        nend = 100
        for _ in range(trial):
            q = np.random.rand(3, 1) * 0.05
            #print(q)
            for _ in range(nend):
                e = np.linalg.norm(xd - self.calc_X(q, xi=1))
                #print(e)
                if e < 1e-3:
                    return q
                else:
                    H = self.hessian(xd, q, xi=1)
                    grad = self.gradient(xd, q, xi=1)
                    
                    #print("H = ", H)
                    #print("grad = ", grad)
                    
                    q += -np.linalg.inv(H) @ grad
        
        return np.zeros((3, 1))
    

    
    






if __name__ == "__main__":
    pass