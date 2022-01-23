import numpy as np

import J_0
import J_1
import J_2
import J_3
import J_4



class DifferentialKinematics:
    """微分運動学"""
    
    def J_0(self, q, xi=1):
        return J_0.f(q, xi)

    def J_1(self, q, xi=1):
        return J_1.f(q, xi)
    
    def J_2(self, q, xi=1):
        return J_2.f(q, xi)
    
    def J_3(self, q, xi=1):
        return J_3.f(q, xi)
    
    def J_4(self, q, xi=1):
        return J_4.f(q, xi)


if __name__ == "__main__":
    
    hoge = DifferentialKinematics()
    
    a = hoge.J_4(
        q = np.zeros((3*5, 1))
    )
    print(a)