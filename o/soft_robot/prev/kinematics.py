import numpy as np
from math import sin, cos, sqrt

from math_utils import *


class Base:
    """ベース"""

    # モーダル同時変換行列のパラメータ
    c1 = 837019575
    c2 = 4133430
    c3 = 32805
    c4 = 486
    c5 = 18
    c6 = 55801305
    c7 = 688905
    c8 = 3645
    c9 = 81
    
    c10 = 279006525
    c11 = 1377810
    c12 = 10935
    c13 = 162

    c14 = 243
    c15 = 2066715

    r = 0.0125
    L0 = 0.15
    
    sq3 = sqrt(3)
    
    
    
    def calc_P(self, q, xi):
        """線形化されたアクチュエータ空間からタスク空間への写像
        
        順運動学
        """
        
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        
        A1 = l1**2 + l2**2 + l3**2 - l1*l2 - l1*l3 - l2*l3
        A2 = 2*l1 - l2 - l3
        A3 = l2 - l3
        A4 = 3*self.L0 + l1 + l2 + l3
        
        x = -(A2 * A1**4 * A4 * xi**10) / ((self.c1 * self.r**9)) + \
            (A2 * A1**3 * A4 * xi**8) / (self.c2 * self.r**7) - \
                (A2 * A1**2 * A4 * xi**6) / (self.c3 * self.r**5) + \
                    (A2 * A1 * A4 * xi**4) / (self.c4 * self.r**3) - \
                        (A2 * A4 * xi**2) / (self.c5 * self.r)
        
        y = -(self.sq3 * A4 * A3 * A1**4 * xi**10) / (self.c1 * self.r**9) + \
            (self.sq3 * A4 * A3 * A1**3 * xi**8) / (self.c2 * self.r**7) - \
                (self.sq3 * A4 * A3 * A1**2 * xi**6) / (self.c3 * self.r**5) + \
                    (self.sq3 * A4 * A1 * A2 * xi**4) / (self.c4 * self.r**3) - \
                        (self.sq3 * A4 * A3 * xi**2) / (self.c5 * self.r)
        
        z = (2 * A1**4 * A4 * xi**9) / (self.c6 * self.r**8) - \
            (4 * A1**3 * A4 * xi**7) / (self.c7 * self.r**6) + \
                (2 * A1**2 * A4 * xi**5) / (self.c8 * self.r**4) - \
                    (2 * A1 *A4 * xi**3) / (self.c9 * self.r**2) + \
                        (A4 * xi) / 3

        return np.array([[x, y, z]]).T
    


    def calc_R(self, q, xi):
        """線形化された回転行列"""
        
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        
        A1 = l1**2 + l2**2 + l3**2 - l1*l2 - l1*l3 - l2*l3
        A2 = 2*l1 - l2 - l3
        A3 = l2 - l3
        A4 = 3*self.L0 + l1 + l2 + l3
        
        R11 = 1 - (A2**2 * A1**4 * xi**10) / (self.c1 * self.r**10) + \
            (A2**2 * A1**3 * xi**8) / (self.c1 * self.r**8) - \
                (A2**2 * A1**2 * xi**6) / (self.c3 * self.r**6) + \
                    (A1 * A2**2 * xi**4) / (self.c4 * self.r**4) - \
                        (A2**2 * xi**2) / (self.c5 * self.r**2)
        
        R12 = (self.sq3 * A2 * A3 * A1**4 * xi**10) / (self.c1 * self.r**10) + \
            (self.sq3 * A2 * A3 * A1**3 * xi**8) / (self.c2 * self.r**8) - \
                (self.sq3 * A2 * A3 * A1**2 * xi**6) / (self.c3 * self.r**6) + \
                    (self.sq3 * A2 * A3 * A1 * xi**4) / (self.c4 * self.r**4) - \
                        (self.sq3 * A2 * A3 * xi**2) / (self.c5 * self.r**2)
        
        R13 = -(2 * A2 * A1**4 * xi**9) / (self.c6 * self.r**9) + \
            (4 * A2 * A1**3 * xi**7) / (self.c7 * self.r**7) - \
                (2 * A2 * A1**2 * xi**5) / (self.c8 * self.r**5) + \
                    (2 * A2 * A1 * xi**3) / (self.c9 * self.r**3) - \
                        (A2 * xi) / (3 * self.r)
        
        R22 = 1 - (A3**2 * A1**4 * xi**10) / (self.c10 * self.r**10) + \
            (A3**2 * A1**3 * xi**8) / (self.c11 * self.r**8) - \
                (A3**2 * A1**2 * xi**6) / (self.c12 * self.r**6) + \
                    (A3**2 * A1 * xi**4) / (self.c13 * self.r**4) - \
                        (A3**2 * xi**2) / (6 * self.r**2)
        
        R23 = -(2*self.sq3 * A3 * A1**4 * xi**9) / (self.c6 * self.r**9) + \
            (4*self.sq3 * A3 * A1**3 * xi**7) / (self.c7 * self.r**7) - \
                (2*self.sq3 * A3 * A1**2 * xi**5) / (self.c8 * self.r**5) + \
                    (2*self.sq3 * A3 * A1 * xi**3) / (self.c9 * self.r**3) - \
                        (self.sq3 * A3 * xi) / (3 * self.r)
        
        R33 = 1 - (2 * xi**2 * A1) / (9 * self.r**2) + \
            (2 * xi**4 * A1**2) / (self.c14 * self.r**4) - \
                (4 * xi**6 * A1**3) / (self.c3 * self.r**6) + \
                    (2 * xi**8 * A1**4) / (self.c15 * self.r**8) - \
                        (4 * xi**10 * A1**5) / (self.c1 * self.r**10)
        
        R21 = R12
        R31 = -R13
        R32 = -R23
        
        return np.array([
            [R11, R12, R13],
            [R21, R22, R23],
            [R31, R32, R33],
        ])


    def calc_MHTM(self, q, xi):
        """モーダル同時変換行列
        
        線形化されたHomogeneous Transformation Matrix
        """
        return np.block([
            [self.calc_R(q, xi), self.calc_P(q, xi)],
            [np.zeros((1, 3)), np.eye(1)],
        ])

    def calc_dPdl1(self, q, xi):
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        
        return np.array([[
            -xi**2*(2*l1 - l2 - l3)/(self.c5*self.r) - 2*xi**2*(3*self.L0 + l1 + l2 + l3)/(self.c5*self.r) + xi**4*(2*l1 - l2 - l3)**2*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) + 2*xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - xi**6*(2*l1 - l2 - l3)*(4*l1 - 2*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - xi**6*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) - 2*xi**6*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + xi**8*(2*l1 - l2 - l3)*(6*l1 - 3*l2 - 3*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + xi**8*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) + 2*xi**8*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - xi**10*(2*l1 - l2 - l3)*(8*l1 - 4*l2 - 4*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - xi**10*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9) - 2*xi**10*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9),
            -self.sq3*xi**2*(l2 - l3)/(self.c5*self.r) + self.sq3*xi**4*(2*l1 - l2 - l3)**2*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + self.sq3*xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) + 2*self.sq3*xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - self.sq3*xi**6*(l2 - l3)*(4*l1 - 2*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - self.sq3*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + self.sq3*xi**8*(l2 - l3)*(6*l1 - 3*l2 - 3*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + self.sq3*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - self.sq3*xi**10*(l2 - l3)*(8*l1 - 4*l2 - 4*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - self.sq3*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9),
            xi/3 - xi**3*(4*l1 - 2*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)/(self.c9*self.r**2) - xi**3*(2*l1**2 - 2*l1*l2 - 2*l1*l3 + 4*l2**2 - 2*l2*l3)/(self.c9*self.r**2) + 2*xi**5*(4*l1 - 2*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c8*self.r**4) + 2*xi**5*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c8*self.r**4) - 4*xi**7*(6*l1 - 3*l2 - 3*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c7*self.r**6) - 4*xi**7*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c7*self.r**6) + 2*xi**9*(8*l1 - 4*l2 - 4*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c6*self.r**8) + 2*xi**9*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c6*self.r**8),
        ]]).T

    def calc_dPdl2(self, q, xi):
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        
        return np.array([[
            -xi**2*(2*l1 - l2 - l3)/(self.c5*self.r) + xi**2*(3*self.L0 + l1 + l2 + l3)/(self.c5*self.r) + xi**4*(-l1 + 4*l2 - l3)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - xi**6*(-2*l1 + 8*l2 - 2*l3)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - xi**6*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + xi**6*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + xi**8*(-3*l1 + 12*l2 - 3*l3)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + xi**8*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - xi**8*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - xi**10*(-4*l1 + 16*l2 - 4*l3)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - xi**10*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9) + xi**10*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9),
            -self.sq3*xi**2*(l2 - l3)/(self.c5*self.r) - self.sq3*xi**2*(3*self.L0 + l1 + l2 + l3)/(self.c5*self.r) + self.sq3*xi**4*(-l1 + 4*l2 - l3)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + self.sq3*xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - self.sq3*xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - self.sq3*xi**6*(l2 - l3)*(-2*l1 + 8*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - self.sq3*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) - self.sq3*xi**6*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + self.sq3*xi**8*(l2 - l3)*(-3*l1 + 12*l2 - 3*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + self.sq3*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) + self.sq3*xi**8*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - self.sq3*xi**10*(l2 - l3)*(-4*l1 + 16*l2 - 4*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - self.sq3*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9) - self.sq3*xi**10*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9),
            xi/3 - xi**3*(-2*l1 + 8*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)/(self.c9*self.r**2) - xi**3*(2*l1**2 - 2*l1*l2 - 2*l1*l3 + 4*l2**2 - 2*l2*l3)/(self.c9*self.r**2) + 2*xi**5*(-2*l1 + 8*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c8*self.r**4) + 2*xi**5*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c8*self.r**4) - 4*xi**7*(-3*l1 + 12*l2 - 3*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c7*self.r**6) - 4*xi**7*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c7*self.r**6) + 2*xi**9*(-4*l1 + 16*l2 - 4*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c6*self.r**8) + 2*xi**9*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c6*self.r**8),
        ]]).T

    def calc_dPdl3(self, q, xi):
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]

        return np.array([[
            -xi**2*(2*l1 - l2 - l3)/(self.c5*self.r) + xi**2*(3*self.L0 + l1 + l2 + l3)/(self.c5*self.r) + xi**4*(-l1 - l2)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - xi**6*(-2*l1 - 2*l2)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - xi**6*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + xi**6*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + xi**8*(-3*l1 - 3*l2)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + xi**8*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - xi**8*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - xi**10*(-4*l1 - 4*l2)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - xi**10*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9) + xi**10*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9),
            -self.sq3*xi**2*(l2 - l3)/(self.c5*self.r) + self.sq3*xi**2*(3*self.L0 + l1 + l2 + l3)/(self.c5*self.r) + self.sq3*xi**4*(-l1 - l2)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + self.sq3*xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - self.sq3*xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - self.sq3*xi**6*(-2*l1 - 2*l2)*(l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - self.sq3*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + self.sq3*xi**6*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + self.sq3*xi**8*(-3*l1 - 3*l2)*(l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + self.sq3*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - self.sq3*xi**8*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - self.sq3*xi**10*(-4*l1 - 4*l2)*(l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - self.sq3*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9) + self.sq3*xi**10*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9),
            xi/3 - xi**3*(-2*l1 - 2*l2)*(3*self.L0 + l1 + l2 + l3)/(self.c9*self.r**2) - xi**3*(2*l1**2 - 2*l1*l2 - 2*l1*l3 + 4*l2**2 - 2*l2*l3)/(self.c9*self.r**2) + 2*xi**5*(-2*l1 - 2*l2)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c8*self.r**4) + 2*xi**5*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c8*self.r**4) - 4*xi**7*(-3*l1 - 3*l2)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c7*self.r**6) - 4*xi**7*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c7*self.r**6) + 2*xi**9*(-4*l1 - 4*l2)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c6*self.r**8) + 2*xi**9*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c6*self.r**8)
        ]]).T


    def calc_dRdl1(self, q, xi):
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        return np.array([
            [
                -xi**2*(8*l1 - 4*l2 - 4*l3)/(self.c5*self.r**2) + xi**4*(2*l1 - l2 - l3)**3/(self.c4*self.r**4) + xi**4*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) - xi**6*(2*l1 - l2 - l3)**2*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c3*self.r**6) - xi**6*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + xi**8*(2*l1 - l2 - l3)**2*(6*l1 - 3*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c1*self.r**8) + xi**8*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**8) - xi**10*(2*l1 - l2 - l3)**2*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**10) - xi**10*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10), 
                -2*self.sq3*xi**2*(l2 - l3)/(self.c5*self.r**2) + self.sq3*xi**4*(l2 - l3)*(2*l1 - l2 - l3)**2/(self.c4*self.r**4) + 2*self.sq3*xi**4*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) - self.sq3*xi**6*(l2 - l3)*(2*l1 - l2 - l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c3*self.r**6) - 2*self.sq3*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + self.sq3*xi**8*(l2 - l3)*(2*l1 - l2 - l3)*(6*l1 - 3*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c2*self.r**8) + 2*self.sq3*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c2*self.r**8) + self.sq3*xi**10*(l2 - l3)*(2*l1 - l2 - l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**10) + 2*self.sq3*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10),
                -2*xi/(3*self.r) + xi**3*(2*l1 - l2 - l3)*(4*l1 - 2*l2 - 2*l3)/(self.c9*self.r**3) + 4*xi**3*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c9*self.r**3) - xi**5*(4*l1 - 2*l2 - 2*l3)**2*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) - 4*xi**5*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c8*self.r**5) + xi**7*(6*l1 - 3*l2 - 3*l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) + 8*xi**7*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c7*self.r**7) - xi**9*(4*l1 - 2*l2 - 2*l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9) - 4*xi**9*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c6*self.r**9)
            ], 
            [
                -2*self.sq3*xi**2*(l2 - l3)/(self.c5*self.r**2) + self.sq3*xi**4*(l2 - l3)*(2*l1 - l2 - l3)**2/(self.c4*self.r**4) + 2*self.sq3*xi**4*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) - self.sq3*xi**6*(l2 - l3)*(2*l1 - l2 - l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c3*self.r**6) - 2*self.sq3*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + self.sq3*xi**8*(l2 - l3)*(2*l1 - l2 - l3)*(6*l1 - 3*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c2*self.r**8) + 2*self.sq3*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c2*self.r**8) + self.sq3*xi**10*(l2 - l3)*(2*l1 - l2 - l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**10) + 2*self.sq3*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10), 
                xi**4*(l2 - l3)**2*(2*l1 - l2 - l3)/(self.c13*self.r**4) - xi**6*(l2 - l3)**2*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c12*self.r**6) + xi**8*(l2 - l3)**2*(6*l1 - 3*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c11*self.r**8) - xi**10*(l2 - l3)**2*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c10*self.r**10),
                2*self.sq3*xi**3*(l2 - l3)*(2*l1 - l2 - l3)/(self.c9*self.r**3) - 2*self.sq3*xi**5*(l2 - l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) + 4*self.sq3*xi**7*(l2 - l3)*(6*l1 - 3*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) - 2*self.sq3*xi**9*(l2 - l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9)
            ],
            [
                2*xi/(3*self.r) - xi**3*(2*l1 - l2 - l3)*(4*l1 - 2*l2 - 2*l3)/(self.c9*self.r**3) - 4*xi**3*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c9*self.r**3) + xi**5*(4*l1 - 2*l2 - 2*l3)**2*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) + 4*xi**5*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c8*self.r**5) - xi**7*(6*l1 - 3*l2 - 3*l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) - 8*xi**7*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c7*self.r**7) + xi**9*(4*l1 - 2*l2 - 2*l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9) + 4*xi**9*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c6*self.r**9),
                -2*self.sq3*xi**3*(l2 - l3)*(2*l1 - l2 - l3)/(self.c9*self.r**3) + 2*self.sq3*xi**5*(l2 - l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) - 4*self.sq3*xi**7*(l2 - l3)*(6*l1 - 3*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) + 2*self.sq3*xi**9*(l2 - l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9),
                -2*xi**2*(2*l1 - l2 - l3)/(9*self.r**2) - 4*xi**6*(6*l1 - 3*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + 2*xi**8*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c15*self.r**8) + 2*xi**4*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c14*self.r**4) - 4*xi**10*(10*l1 - 5*l2 - 5*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10)
            ]
        ])
    
    
    def calc_dRdl2(self, q, xi):
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        return np.array([
            [
                -xi**2*(-4*l1 + 2*l2 + 2*l3)/(self.c5*self.r**2) + xi**4*(-4*l1 + 2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) + xi**4*(-l1 + 2*l2 - l3)*(2*l1 - l2 - l3)**2/(self.c4*self.r**4) - xi**6*(-4*l1 + 2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) - xi**6*(-2*l1 + 4*l2 - 2*l3)*(2*l1 - l2 - l3)**2*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c3*self.r**6) + xi**8*(-4*l1 + 2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**8) + xi**8*(-3*l1 + 6*l2 - 3*l3)*(2*l1 - l2 - l3)**2*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c1*self.r**8) - xi**10*(-4*l1 + 2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10) - xi**10*(-4*l1 + 8*l2 - 4*l3)*(2*l1 - l2 - l3)**2*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**10), 
                self.sq3*xi**2*(l2 - l3)/(self.c5*self.r**2) - self.sq3*xi**2*(2*l1 - l2 - l3)/(self.c5*self.r**2) + self.sq3*xi**4*(l2 - l3)*(-l1 + 2*l2 - l3)*(2*l1 - l2 - l3)/(self.c4*self.r**4) - self.sq3*xi**4*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) + self.sq3*xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) - self.sq3*xi**6*(l2 - l3)*(-2*l1 + 4*l2 - 2*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c3*self.r**6) + self.sq3*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) - self.sq3*xi**6*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + self.sq3*xi**8*(l2 - l3)*(-3*l1 + 6*l2 - 3*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c2*self.r**8) - self.sq3*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c2*self.r**8) + self.sq3*xi**8*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c2*self.r**8) + self.sq3*xi**10*(l2 - l3)*(-4*l1 + 8*l2 - 4*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**10) - self.sq3*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10) + self.sq3*xi**10*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10),
                xi/(3*self.r) + xi**3*(-l1 + 2*l2 - l3)*(4*l1 - 2*l2 - 2*l3)/(self.c9*self.r**3) - 2*xi**3*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c9*self.r**3) - xi**5*(-2*l1 + 4*l2 - 2*l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) + 2*xi**5*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c8*self.r**5) + xi**7*(-3*l1 + 6*l2 - 3*l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) - 4*xi**7*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c7*self.r**7) - xi**9*(-4*l1 + 8*l2 - 4*l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9) + 2*xi**9*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c6*self.r**9),
            ],
            [
                self.sq3*xi**2*(l2 - l3)/(self.c5*self.r**2) - self.sq3*xi**2*(2*l1 - l2 - l3)/(self.c5*self.r**2) + self.sq3*xi**4*(l2 - l3)*(-l1 + 2*l2 - l3)*(2*l1 - l2 - l3)/(self.c4*self.r**4) - self.sq3*xi**4*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) + self.sq3*xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) - self.sq3*xi**6*(l2 - l3)*(-2*l1 + 4*l2 - 2*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c3*self.r**6) + self.sq3*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) - self.sq3*xi**6*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + self.sq3*xi**8*(l2 - l3)*(-3*l1 + 6*l2 - 3*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c2*self.r**8) - self.sq3*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c2*self.r**8) + self.sq3*xi**8*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c2*self.r**8) + self.sq3*xi**10*(l2 - l3)*(-4*l1 + 8*l2 - 4*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**10) - self.sq3*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10) + self.sq3*xi**10*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10),
                -xi**2*(2*l2 - 2*l3)/(6*self.r**2) + xi**4*(l2 - l3)**2*(-l1 + 2*l2 - l3)/(self.c13*self.r**4) + xi**4*(2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c13*self.r**4) - xi**6*(l2 - l3)**2*(-2*l1 + 4*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c12*self.r**6) - xi**6*(2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c12*self.r**6) + xi**8*(l2 - l3)**2*(-3*l1 + 6*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c11*self.r**8) + xi**8*(2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c11*self.r**8) - xi**10*(l2 - l3)**2*(-4*l1 + 8*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c10*self.r**10) - xi**10*(2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c10*self.r**10),
                -self.sq3*xi/(3*self.r) + 2*self.sq3*xi**3*(l2 - l3)*(-l1 + 2*l2 - l3)/(self.c9*self.r**3) + 2*self.sq3*xi**3*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c9*self.r**3) - 2*self.sq3*xi**5*(l2 - l3)*(-2*l1 + 4*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) - 2*self.sq3*xi**5*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c8*self.r**5) + 4*self.sq3*xi**7*(l2 - l3)*(-3*l1 + 6*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) + 4*self.sq3*xi**7*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c7*self.r**7) - 2*self.sq3*xi**9*(l2 - l3)*(-4*l1 + 8*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9) - 2*self.sq3*xi**9*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c6*self.r**9),
            ], 
            [
                -xi/(3*self.r) - xi**3*(-l1 + 2*l2 - l3)*(4*l1 - 2*l2 - 2*l3)/(self.c9*self.r**3) + 2*xi**3*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c9*self.r**3) + xi**5*(-2*l1 + 4*l2 - 2*l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) - 2*xi**5*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c8*self.r**5) - xi**7*(-3*l1 + 6*l2 - 3*l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) + 4*xi**7*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c7*self.r**7) + xi**9*(-4*l1 + 8*l2 - 4*l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9) - 2*xi**9*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c6*self.r**9),
                self.sq3*xi/(3*self.r) - 2*self.sq3*xi**3*(l2 - l3)*(-l1 + 2*l2 - l3)/(self.c9*self.r**3) - 2*self.sq3*xi**3*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c9*self.r**3) + 2*self.sq3*xi**5*(l2 - l3)*(-2*l1 + 4*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) + 2*self.sq3*xi**5*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c8*self.r**5) - 4*self.sq3*xi**7*(l2 - l3)*(-3*l1 + 6*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) - 4*self.sq3*xi**7*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c7*self.r**7) + 2*self.sq3*xi**9*(l2 - l3)*(-4*l1 + 8*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9) + 2*self.sq3*xi**9*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c6*self.r**9),
                -2*xi**2*(-l1 + 2*l2 - l3)/(9*self.r**2) - 4*xi**6*(-3*l1 + 6*l2 - 3*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + 2*xi**8*(-4*l1 + 8*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c15*self.r**8) + 2*xi**4*(-2*l1 + 4*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c14*self.r**4) - 4*xi**10*(-5*l1 + 10*l2 - 5*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10)
            ]
        ])

    def calc_dRdl3(self, q, xi):
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        return np.array([
            [
                -xi**2*(-4*l1 + 2*l2 + 2*l3)/(self.c5*self.r**2) + xi**4*(-4*l1 + 2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) + xi**4*(-l1 - l2 + 2*l3)*(2*l1 - l2 - l3)**2/(self.c4*self.r**4) - xi**6*(-4*l1 + 2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) - xi**6*(-2*l1 - 2*l2 + 4*l3)*(2*l1 - l2 - l3)**2*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c3*self.r**6) + xi**8*(-4*l1 + 2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**8) + xi**8*(-3*l1 - 3*l2 + 6*l3)*(2*l1 - l2 - l3)**2*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c1*self.r**8) - xi**10*(-4*l1 - 4*l2 + 8*l3)*(2*l1 - l2 - l3)**2*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**10) - xi**10*(-4*l1 + 2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10),
                self.sq3*xi**2*(l2 - l3)/(self.c5*self.r**2) + self.sq3*xi**2*(2*l1 - l2 - l3)/(self.c5*self.r**2) + self.sq3*xi**4*(l2 - l3)*(-l1 - l2 + 2*l3)*(2*l1 - l2 - l3)/(self.c4*self.r**4) - self.sq3*xi**4*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) - self.sq3*xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) - self.sq3*xi**6*(l2 - l3)*(-2*l1 - 2*l2 + 4*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c3*self.r**6) + self.sq3*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + self.sq3*xi**6*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + self.sq3*xi**8*(l2 - l3)*(-3*l1 - 3*l2 + 6*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c2*self.r**8) - self.sq3*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c2*self.r**8) - self.sq3*xi**8*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c2*self.r**8) + self.sq3*xi**10*(l2 - l3)*(-4*l1 - 4*l2 + 8*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**10) - self.sq3*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10) - self.sq3*xi**10*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10),
                xi/(3*self.r) + xi**3*(-l1 - l2 + 2*l3)*(4*l1 - 2*l2 - 2*l3)/(self.c9*self.r**3) - 2*xi**3*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c9*self.r**3) - xi**5*(-2*l1 - 2*l2 + 4*l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) + 2*xi**5*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c8*self.r**5) + xi**7*(-3*l1 - 3*l2 + 6*l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) - 4*xi**7*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c7*self.r**7) - xi**9*(-4*l1 - 4*l2 + 8*l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9) + 2*xi**9*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c6*self.r**9)
            ],
            [
                self.sq3*xi**2*(l2 - l3)/(self.c5*self.r**2) + self.sq3*xi**2*(2*l1 - l2 - l3)/(self.c5*self.r**2) + self.sq3*xi**4*(l2 - l3)*(-l1 - l2 + 2*l3)*(2*l1 - l2 - l3)/(self.c4*self.r**4) - self.sq3*xi**4*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) - self.sq3*xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c4*self.r**4) - self.sq3*xi**6*(l2 - l3)*(-2*l1 - 2*l2 + 4*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c3*self.r**6) + self.sq3*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + self.sq3*xi**6*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + self.sq3*xi**8*(l2 - l3)*(-3*l1 - 3*l2 + 6*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c2*self.r**8) - self.sq3*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c2*self.r**8) - self.sq3*xi**8*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c2*self.r**8) + self.sq3*xi**10*(l2 - l3)*(-4*l1 - 4*l2 + 8*l3)*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c1*self.r**10) - self.sq3*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10) - self.sq3*xi**10*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10),
                -xi**2*(-2*l2 + 2*l3)/(6*self.r**2) + xi**4*(-2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c13*self.r**4) + xi**4*(l2 - l3)**2*(-l1 - l2 + 2*l3)/(self.c13*self.r**4) - xi**6*(-2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c12*self.r**6) - xi**6*(l2 - l3)**2*(-2*l1 - 2*l2 + 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c12*self.r**6) + xi**8*(-2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c11*self.r**8) + xi**8*(l2 - l3)**2*(-3*l1 - 3*l2 + 6*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c11*self.r**8) - xi**10*(-2*l2 + 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c10*self.r**10) - xi**10*(l2 - l3)**2*(-4*l1 - 4*l2 + 8*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c10*self.r**10),
                self.sq3*xi/(3*self.r) + 2*self.sq3*xi**3*(l2 - l3)*(-l1 - l2 + 2*l3)/(self.c9*self.r**3) - 2*self.sq3*xi**3*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c9*self.r**3) - 2*self.sq3*xi**5*(l2 - l3)*(-2*l1 - 2*l2 + 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) + 2*self.sq3*xi**5*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c8*self.r**5) + 4*self.sq3*xi**7*(l2 - l3)*(-3*l1 - 3*l2 + 6*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) - 4*self.sq3*xi**7*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c7*self.r**7) - 2*self.sq3*xi**9*(l2 - l3)*(-4*l1 - 4*l2 + 8*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9) + 2*self.sq3*xi**9*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c6*self.r**9)], 
            [
                -xi/(3*self.r) - xi**3*(-l1 - l2 + 2*l3)*(4*l1 - 2*l2 - 2*l3)/(self.c9*self.r**3) + 2*xi**3*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c9*self.r**3) + xi**5*(-2*l1 - 2*l2 + 4*l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) - 2*xi**5*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c8*self.r**5) - xi**7*(-3*l1 - 3*l2 + 6*l3)*(8*l1 - 4*l2 - 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) + 4*xi**7*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c7*self.r**7) + xi**9*(-4*l1 - 4*l2 + 8*l3)*(4*l1 - 2*l2 - 2*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9) - 2*xi**9*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c6*self.r**9), -self.sq3*xi/(3*self.r) - 2*self.sq3*xi**3*(l2 - l3)*(-l1 - l2 + 2*l3)/(self.c9*self.r**3) + 2*self.sq3*xi**3*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c9*self.r**3) + 2*self.sq3*xi**5*(l2 - l3)*(-2*l1 - 2*l2 + 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c8*self.r**5) - 2*self.sq3*xi**5*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c8*self.r**5) - 4*self.sq3*xi**7*(l2 - l3)*(-3*l1 - 3*l2 + 6*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c7*self.r**7) + 4*self.sq3*xi**7*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c7*self.r**7) + 2*self.sq3*xi**9*(l2 - l3)*(-4*l1 - 4*l2 + 8*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c6*self.r**9) - 2*self.sq3*xi**9*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c6*self.r**9),
                -2*xi**2*(-l1 - l2 + 2*l3)/(9*self.r**2) - 4*xi**6*(-3*l1 - 3*l2 + 6*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**2/(self.c3*self.r**6) + 2*xi**8*(-4*l1 - 4*l2 + 8*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**3/(self.c15*self.r**8) + 2*xi**4*(-2*l1 - 2*l2 + 4*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)/(self.c14*self.r**4) - 4*xi**10*(-5*l1 - 5*l2 + 10*l3)*(l1**2 - l1*l2 - l1*l3 + l2**2 - l2*l3 + l3**2)**4/(self.c1*self.r**10)
            ]
        ])


class OneSection(Base):
    """1つセクションのローカル変数"""
    
    def __init__(self, n, n_step):
        
        self.n = n  # セクション番号
        
        self.J_OMEGAs = None
        self.J_omegas = None
        self.J_vs = None
    
    
    def P()




class AllSection:
    """アーム全体の運動学"""
    
    def __init__(self, N):
        
        self.N = N  # セクションの数
        self.n_step = 10
        self.set_section()
    
    
    def set_section(self,):
        """ローカルセクションを追加"""
        self.sections = [
            OneSection(i, self.n_step) for i in range(self.N)
        ]
    
    
    def update_all(self, q_all, q_dot_all):
        """全部更新"""
        self.q_all = q_all  # 縦ベクトル
        self.q_dot_all = q_dot_all
    
    
    def update_local(self,):
        """ローカル位置，回転行列を更新"""
        for i in range(self.N):
            self.sections[i].update_local_state(self.q_all[i:i+3, :])
    
    
    def update_J_OMEGA_ij(self,):
        """J_OMEGAを更新"""
        
        for k in range(self.N):
            print("k = ", k)
            
            J_OMEGA_ijs_all = []
            for l in range(self.n_step):  # xiの一個一個の分を順番に計算
                J_OMEGA_ijs = []
                for i in range(k+1):
                    print("i = ", i)
                    Ri = self.sections[i].Rs[l]
                    J_OMEGA_ij = []
                    for j in range(3):
                        
                        if i == k:
                            print("ketu")
                            dRidlj = self.sections[i].dRdls[l][j]
                            J_OMEGA_ij.append(Ri.T @ dRidlj)
                        else:
                            print("mae")
                            print(self.sections[i].J_OMEGAs)
                            J_OMEGGA_prev = self.sections[i-1].J_OMEGAs[-1][i][j]
                            J_OMEGA_ij.append(Ri.T @ J_OMEGGA_prev @ Ri)
                    
                    J_OMEGA_ijs.append(J_OMEGA_ij)
                J_OMEGA_ijs_all.append(J_OMEGA_ijs)
            
            self.sections[k].J_OMEGAs = J_OMEGA_ijs_all

        
        
        return
    
    
    
    def update_J_v_ij(self,):
        pass


if __name__ == "__main__":
    N = 3
    q_all = np.zeros((3*N, 1))
    q_dot_all = np.zeros((3*N, 1))

    
    hoge = AllSection(N)
    hoge.update_all(q_all, q_dot_all)
    hoge.update_local()
    hoge.update_J_OMEGA_ij()
