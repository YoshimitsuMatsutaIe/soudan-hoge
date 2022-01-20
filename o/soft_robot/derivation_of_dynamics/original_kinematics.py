"""低次多変量テイラー展開"""

import sympy as sy
from sympy import sin, cos, sqrt, atan2, pi


class Local:
    """オリジナルの式"""
    r = 0.0125
    L0 = 0.15
    sq3 = sqrt(3)
    
    def mapping_from_actuator_to_configuration(self, q):
        """アクチュエータ空間から配置空間への写像"""
        
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        
        A1 = l1**2 + l2**2 + l3**2 - \
            l1*l2 - l1*l3 - l2*l3
        A2 = 2*l1 - l2 - l3
        A3 = l2 - l3
        A4 = 3*self.L0 + l1 + l2 + l3
        
        #print("A1 =", A1)
        # if A1 < 1e-6:
        #     lam = A4 * self.r / 2*1e-6
        # else:
        #     lam = A4 * self.r / 2*sqrt(A1)
        lam = A4 * self.r / 2*sqrt(A1)
        phi = 2*sqrt(A1) / 3*self.r
        theta = atan2(self.sq3 * (-A3) / (-A2), 1)
        
        # if phi <= 0 or phi > 2*pi:
        #     print("phiが範囲外!")
        
        # if theta <= -pi or theta >= pi:
        #     print("thetaが範囲外")
        
        return sy.Matrix([[lam, phi, theta]]).T


    def mapping_from_configration_to_task_p(self, c, xi):
        """配置空間からタスク空間pへの写像"""
        
        lam = c[0, 0]
        phi = c[1, 0]
        theta = c[2, 0]
        
        return sy.Matrix([
            [lam * cos(theta) * (1 - cos(xi * phi))],
            [lam * sin(theta) * (1 - cos(xi * phi))],
            [lam * sin(xi * phi)],
        ])


    def mapping_from_configration_to_task_R(self, c, xi):
        """配置空間からタスク空間Rへの写像"""
        
        #lam = c[0, 0]
        phi = c[1, 0]
        theta = c[2, 0]
        
        R11 = cos(theta)**2 * cos(xi*phi) + sin(theta)**2
        R12 = sin(theta) * cos(theta) * (cos(xi*phi) - 1)
        R13 = cos(theta) * sin(xi*phi)
        R21 = R12
        R22 = sin(theta)**2 * cos(xi*phi) + cos(theta)**2
        R23 = sin(theta) * sin(xi*phi)
        R31 = -R13
        R32 = -R23
        R33 = cos(xi*phi)
        
        return sy.Matrix([
            [R11, R12, R13],
            [R21, R22, R23],
            [R31, R32, R33],
        ])


    def P(self, q, xi):
        """アクチュエータ空間からタスク空間への写像"""
        
        c = self.mapping_from_actuator_to_configuration(q)
        #print("c = ", c)
        return self.mapping_from_configration_to_task_p(c, xi)



    def R(self, q, xi):
        """アクチュエータ空間からタスク空間への写像"""
        
        c = self.mapping_from_actuator_to_configuration(q)
        return self.mapping_from_configration_to_task_R(c, xi)
