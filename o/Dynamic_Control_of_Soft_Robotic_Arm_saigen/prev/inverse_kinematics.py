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
        """所望のアクチュエータ位置
        
        逆運動学
        """
        
        def f(q):
            q = np.array([q]).T
            z = self.linearized_mapping_from_actuator_to_task_p(q, xi) - pd
            return np.ravel(z).tolist()
    
        ans = fsolve(f, [0, 0, 0])
        return np.array([ans]).T
    
    
    def calc_qd_dot(self, qd, pd_dot, xi=1):
        """所望のアクチュエータ速度"""
        return np.linalg.pinv(self.linearized_jacobian_dpdq(qd, xi)) @ pd_dot
    

    def calc_qd_dot_dot(self, qd, qd_dot, pd_dot_dot, xi=1):
        """所望のアクチュエータ加速度"""
        return np.linalg.pinv(self.linearized_jacobian_dpdq(qd, xi)) @ \
            (pd_dot_dot - self.linearized_jacobian_dpdq_dot(qd, qd_dot, xi) @ qd_dot)
    
    
    def calc_desired_actuator_state(
        self,
        pd_func, pd_dot_func, pd_dot_dot_func,t
    ):
        """時刻tにおける所望のアクチュエータ空間上の状態を計算"""
        qd = self.calc_qd_from_pd_using_fsolve(pd_func(t), xi=1)
        qd_dot = self.calc_qd_dot(qd, pd_dot_func(t), xi=1)
        qd_dot_dot = self.calc_qd_dot_dot(qd, qd_dot, pd_dot_dot_func(t), xi=1)
        
        return qd, qd_dot, qd_dot_dot
    
    
    def calc_desired_actuator_state_all(
        self,
        pd_func, pd_dot_func, pd_dot_dot_func,
        TIME_SPAN, TIME_INTERVAL,
    ):
        """全部計算"""
        
        qd = np.zeros((int(TIME_SPAN / TIME_INTERVAL), 3))
        qd_dot = np.zeros_like(qd)
        qd_dot_dot = np.zeros_like(qd)
        
        for i, t in enumerate(np.arange(0, TIME_SPAN, TIME_INTERVAL)):
            _qd, _qd_dot, _qd_dot_dot = self.calc_desired_actuator_state(
                pd_func, pd_dot_func, pd_dot_dot_func, t
            )
            qd[i:i+1, :] = _qd.T
            qd_dot[i:i+1, :] = _qd_dot.T
            qd_dot_dot[i:i+1, :] = _qd_dot_dot.T
        
        return qd, qd_dot, qd_dot_dot
    
    
    def calc_qd_from_pd_1(self, xd, xi=1,):
        """逆運動学を解く
        
        fsolve使わない実装
        """
        
        
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