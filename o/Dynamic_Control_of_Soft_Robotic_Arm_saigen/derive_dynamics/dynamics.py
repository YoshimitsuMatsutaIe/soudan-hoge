"""ダイナミクスを導出"""

import sympy as sy
from sympy import sqrt


import kinematics

def operate_T2(A):
    """T2演算"""
    return A[0,0] + A[1,1]


def operate_V(A):
    """ベクトル化演算"""
    V_A = A.reshape(len(A), 1)
    return V_A


def operate_tilde(A):
    """行列の最後の列を0にする"""
    m, _ = A.shape
    tilde_A = A
    for i in range(m):
        tilde_A[i, -1] = 0
    return tilde_A


class Dynamics(kinematics.Global):
    """動力学の導出
    
    全セクションで同一のパラメータであると仮定  
    """
    
    m = 1.0
    Ixx = 1.0
    g = sy.Matrix([[0, 0, -9.81]]).T
    
    
    def __init__(self, N):
        
        super().__init__(N)
        self.set_M_omega()





    def set_M_omega(self,):
        """回転方向の慣性行列をセット"""
        
        
        def M_omega(i, j, k):
            
            if j < i and k < i:
                Mijk = self.Ixx * \
                    operate_T2(
                        self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T *\
                            self.J_OMEGA_s[i-1][:, 3*k:3*k+3].T
                    )
                return Mijk.subs(self.xi_large[i, 0], 1)
            
            elif j < i and k == i:
                
        
        
        M_omega_s = []
        for i in range(self.N):
            print("i = ", i)
            M_omega_s_i = []
            
            for j in range(self.N):
                print("j = ", j)
                M_omega_s_ij = []
                
                for k in range(self.N):
                    print("k = ", k)
                    
                    M_omega_s_ij.append(M_omega(i, j, k))
                M_omega_s_i.append(M_omega_s_ij)
            M_omega_s.append(sy.Matrix(M_omega_s_i))
        
        
        self.M_omega_s = M_omega_s






if __name__ == "__main__":
    
    
    hoge = Dynamics(3)