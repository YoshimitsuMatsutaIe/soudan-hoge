"""ダイナミクスを導出"""

import sympy as sy
from sympy import sqrt


import kinematics




class Dynamics(kinematics.Global):
    
    
    def __init__(self, N):
        
        super().__init__(N)
        self.set_M_omega()


    def set_M_omega(self,):
        """回転方向の慣性行列をセット"""
        
        
        print(Q)
        pass






if __name__ == "__main__":
    
    
    hoge = Dynamics(3)