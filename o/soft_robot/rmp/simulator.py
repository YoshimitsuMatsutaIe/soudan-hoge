import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as anm

from kinematics import Kinematics
from differential_kinematics import DifferentialKinematics




class Simulator:
    
    TIME_SPAN = 10  # デフォルト
    TIME_INTERVAL = 0.01  # デフォルト
    
    
    def __init__(self, N=3):
        
        self.N = N  # セクションの数
        
        self.kim = Kinematics()
        self.diff_kim = DifferentialKinematics()
        
        pass
    
    
    def X_dot(t, X):
        """scipyで使う微分方程式"""
        
        
        
    
    
    def run(self, TIME_SPAN=None, TIME_INTERVAL=None):
        """シミュレーション実行"""
        
        if TIME_SPAN is not None:
            self.TIME_SPAN = TIME_SPAN
        
        if TIME_INTERVAL is not None:
            self.TIME_INTERVAL = TIME_INTERVAL
        
        
        sol = integrate.solve_ivp(
            
        )
        
        







if __name__ == "__main__":
    
    hoge = Simulator()
    
    