import numpy as np
from math import sin, cos, sqrt


def sknew_symmetric(a):
    """ベクトル -> 歪対称行列"""
    return np.array([
        [0, -a[2,0], a[1,0]],
        [a[2,0], 0, -a[0,0]],
        [-a[1,0], a[0,0], 0],
    ])

def inv_sknew_symmetric(A):
    """歪対称行列 -> ベクトル"""
    return np.array([[A[2,1], A[0,2], A[1,0]]]).T



class Kinematics:
    """ソフトロボットの運動学"""

    def __init__(self,):
        
        # 線形化した後の様々なパラメータ
        self.c1 = 837019575
        self.c2 = 4133430
        self.c3 = 32805
        self.c4 = 486
        self.c5 = 18
        self.c6 = 55801305
        self.c7 = 688905
        self.c8 = 3645
        self.c9 = 81
        
        self.c10 = 279006525
        self.c11 = 1377810
        self.c12 = 10935
        self.c13 = 162

        self.c14 = 243
        self.c15 = 2066715

        self.r = 0.0125
        self.L0 = 0.15


    def mapping_from_actuator_to_configuration(self, q, xi):
        """アクチュエータ空間から配置空間への写像"""
        
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        
        A1 = l1**2 + l2**2 + l3**2 - \
            l1*l2 - l1*l3 - l2*l3
        A2 = 2*l1 - l2 - l3
        A3 = l2 - l3
        A4 = 3*self.L0 + l1 + l2 + l3
        
        lam = A4 / 2*sqrt(A1)
        phi = 2*sqrt(A1) / 3*self.r
        theta = np.arctan(sqrt(3) * (-A3) / (-A2))
        
        if phi <= 0 or phi > 2*np.pi:
            print("phiが範囲外!")
        
        if theta <= -np.pi or theta >= np.pi:
            print("thetaが範囲外")
        
        return np.array([[lam, phi, theta]]).T


    def mapping_from_configration_to_task(self, c, xi):
        """配置空間からタスク空間への写像"""
        
        lam = c[0, 0]
        phi = c[1, 0]
        theta = c[2, 0]
        
        return np.array([
            [lam * cos(theta) * (1 - cos(xi * phi))],
            [lam * sin(theta) * (1 - cos(xi * phi))],
            [lam * sin(xi * phi)],
        ])


    def mapping_from_actuator_to_task(self, q, xi):
        """アクチュエータ空間からタスク空間への写像"""
        
        c = self.mapping_from_actuator_to_configuration(q, xi)
        x = self.mapping_from_configration_to_task(c, xi)
        
        return x


    def calc_linearized_X(self, q, xi):
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
        
        y = -(sqrt(3) * A4 * A3 * A1**4 * xi**10) / (self.c1 * self.r**9) + \
            (sqrt(3) * A4 * A3 * A1**3 * xi**8) / (self.c2 * self.r**7) - \
                (sqrt(3) * A4 * A3 * A1**2 * xi**6) / (self.c3 * self.r**5) + \
                    (sqrt(3) * A4 * A1 * A2 * xi**4) / (self.c4 * self.r**3) - \
                        (sqrt(3) * A4 * A3 * xi**2) / (self.c5 * self.r)
        
        z = (2 * A1**4 * A4 * xi**9) / (self.c6 * self.r**8) - \
            (4 * A1**3 * A4 * xi**7) / (self.c7 * self.r**6) + \
                (2 * A1**2 * A4 * xi**5) / (self.c8 * self.r**4) - \
                    (2 * A1 *A4 * xi**3) / (self.c9 * self.r**2) + \
                        (A4 * xi) / 3

        return np.array([[x, y, z]]).T


    def calc_linearized_R(self, q, xi):
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
        
        R12 = (sqrt(3) * A2 * A3 * A1**4 * xi**10) / (self.c1 * self.r**10) + \
            (sqrt(3) * A2 * A3 * A1**3 * xi**8) / (self.c2 * self.r**8) - \
                (sqrt(3) * A2 * A3 * A1**2 * xi**6) / (self.c3 * self.r**6) + \
                    (sqrt(3) * A2 * A3 * A1 * xi**4) / (self.c4 * self.r**4) - \
                        (sqrt(3) * A2 * A3 * xi**2) / (self.c5 * self.r**2)
        
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
        
        R23 = -(2*sqrt(3) * A3 * A1**4 * xi**9) / (self.c6 * self.r**9) + \
            (4*sqrt(3) * A3 * A1**3 * xi**7) / (self.c7 * self.r**7) - \
                (2*sqrt(3) * A3 * A1**2 * xi**5) / (self.c8 * self.r**5) + \
                    (2*sqrt(3) * A3 * A1 * xi**3) / (self.c9 * self.r**3) - \
                        (sqrt(3) * A3 * xi) / (3 * self.r)
        
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


    def MHTM(self, q, xi):
        """Modal Homogeneous Transformation Matrix"""
        
        p = self.calc_X(q, xi)
        R = self.calc_R(q, xi)
        
        return np.block([
            [R, p],
            [np.zeros((1, 3)), np.eye(1)],
        ])


    def J_oemga(q, i, j):
        
        
        
        
        
        return 


    def calc_Jacobian(self, q, xi):
        """ヤコビ行列
        
        タスク写像XのアクチュエータベクトルLによる偏微分
        """
        
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        
        return np.array([
            [
                -xi**2*(2*l1 - l2 - l3)/(self.c5*self.r) - 2*xi**2*(3*self.L0 + l1 + l2 + l3)/(self.c5*self.r) + xi**4*(2*l1 - l2 - l3)**2*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) + 2*xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - xi**6*(2*l1 - l2 - l3)*(4*l1 - 2*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - xi**6*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) - 2*xi**6*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + xi**8*(2*l1 - l2 - l3)*(6*l1 - 3*l2 - 3*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + xi**8*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) + 2*xi**8*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - xi**10*(2*l1 - l2 - l3)*(8*l1 - 4*l2 - 4*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - xi**10*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9) - 2*xi**10*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9),
                -xi**2*(2*l1 - l2 - l3)/(self.c5*self.r) + xi**2*(3*self.L0 + l1 + l2 + l3)/(self.c5*self.r) + xi**4*(-l1 + 4*l2 - l3)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - xi**6*(-2*l1 + 8*l2 - 2*l3)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - xi**6*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + xi**6*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + xi**8*(-3*l1 + 12*l2 - 3*l3)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + xi**8*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - xi**8*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - xi**10*(-4*l1 + 16*l2 - 4*l3)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - xi**10*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9) + xi**10*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9),
                -xi**2*(2*l1 - l2 - l3)/(self.c5*self.r) + xi**2*(3*self.L0 + l1 + l2 + l3)/(self.c5*self.r) + xi**4*(-l1 - l2)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - xi**6*(-2*l1 - 2*l2)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - xi**6*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + xi**6*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + xi**8*(-3*l1 - 3*l2)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + xi**8*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - xi**8*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - xi**10*(-4*l1 - 4*l2)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - xi**10*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9) + xi**10*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9)],
            [
                -sqrt(3)*xi**2*(l2 - l3)/(self.c5*self.r) + sqrt(3)*xi**4*(2*l1 - l2 - l3)**2*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + sqrt(3)*xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) + 2*sqrt(3)*xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - sqrt(3)*xi**6*(l2 - l3)*(4*l1 - 2*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - sqrt(3)*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + sqrt(3)*xi**8*(l2 - l3)*(6*l1 - 3*l2 - 3*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + sqrt(3)*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - sqrt(3)*xi**10*(l2 - l3)*(8*l1 - 4*l2 - 4*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - sqrt(3)*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9),
                -sqrt(3)*xi**2*(l2 - l3)/(self.c5*self.r) - sqrt(3)*xi**2*(3*self.L0 + l1 + l2 + l3)/(self.c5*self.r) + sqrt(3)*xi**4*(-l1 + 4*l2 - l3)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + sqrt(3)*xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - sqrt(3)*xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - sqrt(3)*xi**6*(l2 - l3)*(-2*l1 + 8*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - sqrt(3)*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) - sqrt(3)*xi**6*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + sqrt(3)*xi**8*(l2 - l3)*(-3*l1 + 12*l2 - 3*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + sqrt(3)*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) + sqrt(3)*xi**8*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - sqrt(3)*xi**10*(l2 - l3)*(-4*l1 + 16*l2 - 4*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - sqrt(3)*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9) - sqrt(3)*xi**10*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9),
                -sqrt(3)*xi**2*(l2 - l3)/(self.c5*self.r) + sqrt(3)*xi**2*(3*self.L0 + l1 + l2 + l3)/(self.c5*self.r) + sqrt(3)*xi**4*(-l1 - l2)*(2*l1 - l2 - l3)*(3*self.L0 + l1 + l2 + l3)/(self.c4*self.r**3) + sqrt(3)*xi**4*(2*l1 - l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - sqrt(3)*xi**4*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c4*self.r**3) - sqrt(3)*xi**6*(-2*l1 - 2*l2)*(l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c3*self.r**5) - sqrt(3)*xi**6*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + sqrt(3)*xi**6*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c3*self.r**5) + sqrt(3)*xi**8*(-3*l1 - 3*l2)*(l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c2*self.r**7) + sqrt(3)*xi**8*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - sqrt(3)*xi**8*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c2*self.r**7) - sqrt(3)*xi**10*(-4*l1 - 4*l2)*(l2 - l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c1*self.r**9) - sqrt(3)*xi**10*(l2 - l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9) + sqrt(3)*xi**10*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c1*self.r**9)],
            [
                xi/3 - xi**3*(4*l1 - 2*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)/(self.c9*self.r**2) - xi**3*(2*l1**2 - 2*l1*l2 - 2*l1*l3 + 4*l2**2 - 2*l2*l3)/(self.c9*self.r**2) + 2*xi**5*(4*l1 - 2*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c8*self.r**4) + 2*xi**5*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c8*self.r**4) - 4*xi**7*(6*l1 - 3*l2 - 3*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c7*self.r**6) - 4*xi**7*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c7*self.r**6) + 2*xi**9*(8*l1 - 4*l2 - 4*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c6*self.r**8) + 2*xi**9*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c6*self.r**8),
                xi/3 - xi**3*(-2*l1 + 8*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)/(self.c9*self.r**2) - xi**3*(2*l1**2 - 2*l1*l2 - 2*l1*l3 + 4*l2**2 - 2*l2*l3)/(self.c9*self.r**2) + 2*xi**5*(-2*l1 + 8*l2 - 2*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c8*self.r**4) + 2*xi**5*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c8*self.r**4) - 4*xi**7*(-3*l1 + 12*l2 - 3*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c7*self.r**6) - 4*xi**7*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c7*self.r**6) + 2*xi**9*(-4*l1 + 16*l2 - 4*l3)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c6*self.r**8) + 2*xi**9*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c6*self.r**8),
                xi/3 - xi**3*(-2*l1 - 2*l2)*(3*self.L0 + l1 + l2 + l3)/(self.c9*self.r**2) - xi**3*(2*l1**2 - 2*l1*l2 - 2*l1*l3 + 4*l2**2 - 2*l2*l3)/(self.c9*self.r**2) + 2*xi**5*(-2*l1 - 2*l2)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)/(self.c8*self.r**4) + 2*xi**5*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c8*self.r**4) - 4*xi**7*(-3*l1 - 3*l2)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**2/(self.c7*self.r**6) - 4*xi**7*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c7*self.r**6) + 2*xi**9*(-4*l1 - 4*l2)*(3*self.L0 + l1 + l2 + l3)*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**3/(self.c6*self.r**8) + 2*xi**9*(l1**2 - l1*l2 - l1*l3 + 2*l2**2 - l2*l3)**4/(self.c6*self.r**8)
            ]
        ])



if __name__ == "__main__":
    pass