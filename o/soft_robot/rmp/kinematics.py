import numpy as np

import Phi_0
import Phi_1
import Phi_2
import Phi_3
import Phi_4

class Kinematics:
    """運動学"""
    
    
    def Phi_0(self, q, xi=1):
        return Phi_0.f(q, xi)

    def Phi_1(self, q, xi=1):
        return Phi_1.f(q, xi)

    def Phi_2(self, q, xi=1):
        return Phi_2.f(q, xi)
    
    def Phi_3(self, q, xi=1):
        return Phi_3.f(q, xi)
    
    def Phi_4(self, q, xi=1):
        return Phi_4.f(q, xi)





if __name__ == "__main__":
    hoge = Kinematics()
    
    
    a = hoge.Phi_4(np.zeros((15,1)))
    print(a)