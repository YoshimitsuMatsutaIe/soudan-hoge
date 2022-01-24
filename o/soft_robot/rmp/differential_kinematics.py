import numpy as np

import J_0
import J_1
import J_2
import J_3
import J_4



class DifferentialKinematics:
    """微分運動学"""

    def __init__(self, N):
        self.N = N


    def J(self, n, q, xi=1):
        if n == 1:
            return J_0.f(q, xi)[:, :self.N*3]
        elif n == 2:
            return J_1.f(q, xi)[:, :self.N*3]
        elif n == 3:
            return J_2.f(q, xi)[:, :self.N*3]
        elif n == 4:
            return J_3.f(q, xi)[:, :self.N*3]
        elif n == 5:
            return J_4.f(q, xi)[:, :self.N*3]


if __name__ == "__main__":
    
    hoge = DifferentialKinematics(3)
    
    a = hoge.J_4(
        q = np.zeros((3*5, 1))
    )
    print(a)