import numpy as np
import matplotlib.pyplot as plt


import Phi_0
import Phi_1
import Phi_2
import Phi_3
import Phi_4

class Kinematics:
    """運動学"""
    
    xi_large = np.arange(0, 1, 0.01)
    
    def Phi(self, n, q, xi=1):
        if n == 1:
            return Phi_0.f(q, xi)
        elif n == 2:
            return Phi_1.f(q, xi)
        elif n == 3:
            return Phi_2.f(q, xi)
        elif n == 4:
            return Phi_3.f(q, xi)
        elif n == 5:
            return Phi_4.f(q, xi)
        
        return None


class KinematiOriginal:
    """線形化なし運動学
    
    ・基本使わない
    """
    
    def __init__(self,):

        self.r = 0.0125
        self.L0 = 0.15
        
        self.sq3 = np.sqrt(3)


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
        

        lam = A4 * self.r / 2*np.sqrt(A1)
        phi = 2*np.sqrt(A1) / 3*self.r
        theta = np.arctan2(self.sq3 * (-A3) / (-A2), 1)
        
        if phi <= 0 or phi > 2*np.pi:
            print("phiが範囲外!")
        
        if theta <= -np.pi or theta >= np.pi:
            print("thetaが範囲外")
        #print(lam, phi, theta)
        return np.array([[lam, phi, theta]]).T


    def mapping_from_configration_to_task_p(self, c, xi):
        """配置空間からタスク空間pへの写像"""
        
        lam = c[0, 0]
        phi = c[1, 0]
        theta = c[2, 0]
        
        return np.array([
            [lam * np.cos(theta) * (1 - np.cos(xi * phi))],
            [lam * np.sin(theta) * (1 - np.cos(xi * phi))],
            [lam * np.sin(xi * phi)],
        ])



    def mapping_from_actuator_to_task_p(self, q, xi=1):
        """アクチュエータ空間からタスク空間への写像"""
        
        c = self.mapping_from_actuator_to_configuration(q, xi)
        #print("c = ", c)
        x = self.mapping_from_configration_to_task_p(c, xi)
        #print("x = ", x)
        return x



if __name__ == "__main__":
    hoge = Kinematics()
    
    def ikko_dake(l1, l2, l3):
        return hoge.Phi(
                0,
                q = np.array([[
                    l1, l2, l3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                ]]).T
            )
    
    
    l = np.arange(-0.07, 0.07, 0.001)
    
    p = [np.linalg.norm(ikko_dake(l1, 0.001, 0)) for l1 in l]
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(l, p, label="liner_p_norm")
    ax.set_xlabel("l1")
    ax.set_ylabel("[m]")
    
    
    hoge2 = KinematiOriginal()
    
    p = [np.linalg.norm(hoge2.mapping_from_actuator_to_task_p(np.array([[l1, 0.01, 0.6]]).T)) for l1 in l]
    print(p)
    ax.plot(l, p, label="origi_p_norm", linestyle = "dashed")
    
    o = [0.15 for _ in l]
    ax.plot(l, o, label="L0")
    
    ax.legend()
    ax.grid()
    plt.show()