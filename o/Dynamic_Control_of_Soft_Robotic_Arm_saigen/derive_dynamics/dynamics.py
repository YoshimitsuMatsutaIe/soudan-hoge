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
    """行列の最後の列を0埋め"""
    m, _ = A.shape
    tilde_A = A
    for i in range(m):
        tilde_A[i, -1] = 0
    return tilde_A


class Dynamics(kinematics.Global):
    """動力学の導出
    
    全セクションで同一のパラメータであると仮定  
    """
    
    m = 0.13
    Ixx = 1.0
    g = sy.Matrix([[0, 0, -9.81]]).T
    
    
    def __init__(self, N):
        """
        
        Parameters
        ---
        N : int
            セクションの数
        
        N = 3くらいでも計算にかなり時間がかかる
        """
        
        super().__init__(N)
        
        self.set_M_omega()
        self.set_M_omega_dot()





    def set_M_omega(self,):
        """回転方向の慣性行列をセット"""
        print("computing M_omega ...")
        
        def M_omega(i, j, k):
            
            if j < 3*i and k < 3*i:
                Mijk = self.Ixx * \
                    operate_T2(
                        self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T *\
                            self.J_OMEGA_s[i-1][:, 3*k:3*k+3].T
                    )
                return Mijk.subs(self.xi_large[i, 0], 1)
            
            elif j < 3*i and 3*i <= k <= 3*i+2:
                Ri = self.R_s[i]
                Ri_diff_k = sy.diff(Ri, self.q_large[k, 0])
                z = operate_tilde(Ri) * operate_tilde(Ri_diff_k).T
                A = sy.integrate(z, self.xi_large[i, 0])
                A = operate_V(A).T
                
                B = operate_V(self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T)
                B.subs(self.xi_large[i, 0], 1)
                
                return self.Ixx * A * B
            
            elif 3*i <= j <= 3*i+2 and 3*i <= k <= 3*i+2:
                Ri = self.R_s[i]
                Ri_diff_j = sy.diff(Ri, self.q_large[j, 0])
                Ri_diff_k = sy.diff(Ri, self.q_large[k, 0])
                z = Ri_diff_j.T * Ri_diff_k
                A = sy.integrate(z, self.xi_large[i, 0])
                A = operate_T2(A)
                
                return self.Ixx * A
            
            else:
                return 0

        
        M_omega_s = []
        for i in range(self.N):
            #print("i = ", i)
            M_omega_s_i = []
            
            for j in range(3*self.N):
                #print("j = ", j)
                M_omega_s_ij = []
                
                for k in range(3*self.N):
                    #print("k = ", k)
                    
                    M_omega_s_ij.append(M_omega(i, j, k))
                M_omega_s_i.append(M_omega_s_ij)
            M_omega_s.append(sy.Matrix(M_omega_s_i))
        
        
        self.M_omega_s = M_omega_s
        
        print("computing M_omega done!")


    def set_M_omega_dot(self,):
        """回転方向慣性行列の微分をセット"""
        print("computing M_omega_dot ...")
        
        
        def M_omega_dot(i, j, k, s):
            if j < 3*i and k < 3*i:
                if s < 3*i:
                    A = self.H_OMEGA_s[i-1][3*j:3*j+3, 3*s:3*s+3].T *\
                        self.J_OMEGA_s[i-1][:, 3*k:3*k+3]
                    B = self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T *\
                        self.H_OMEGA_s[i-1][3*k:3*k+3, 3*s:3*s+3]
                    Z = self.Ixx * operate_T2(A + B)
                    return Z.subs(self.xi_large[i, 0], 1)
                
                elif 3*i <= s <= 3*i+2:
                    return 0
                
                else:
                    return 0
            
            
            elif j < 3*i and 3*i <= k <= 3*i+2:
                
                Ri = self.R_s[i]
                Ri_diff_k = sy.diff(Ri, self.q_large[k, 0])
                RR_T = operate_tilde(Ri) * operate_tilde(Ri_diff_k).T
                
                if s < 3*i:
                    A = sy.integrate(RR_T, self.xi_large[i, 0])
                    A = operate_V(A).T
                    
                    B = operate_V(self.H_OMEGA_s[i-1][3*j:3*j+3, 3*s:3*s+3].T)
                    
                    Z = self.Ixx * A * B
                    return Z.subs(self.xi_large[i, 0], 1)
                
                elif 3*i <= s <= 3*i+2:
                    RR_T_diff_s = sy.diff(RR_T, self.q_large[s, 0])
                    A = sy.integrate(RR_T_diff_s, self.xi_large[i, 0])
                    A = operate_V(A).T
                    
                    B = operate_V(self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T)
                    
                    Z = self.Ixx * A * B
                    return Z.subs(self.xi_large[i, 0], 1)
                
                else:
                    return 0
            
            elif 3*i <= j <= 3*i+2 and 3*i <= k <= 3*i+2:
                if s < 3*i:
                    return 0
                
                elif 3*i <= s <= 3*i+2:
                    Ri = self.R_s[i]
                    Ri_diff_j = sy.diff(Ri, self.q_large[j, 0])
                    Ri_diff_k = sy.diff(Ri, self.q_large[k, 0])
                    z = Ri_diff_j.T * Ri_diff_k
                    A = sy.integrate(z, self.xi_large[i, 0])
                    A = operate_T2(A)
                    A = sy.diff(A, self.q_large[s, 0])
                    return self.Ixx * A
                
                else:
                    return 0
            
            else:
                return 0

        
        M_omega_dot_s = []
        
        for s in range(3*self.N):
            print("s = ", s)
            M_omega_dot_s_diff_by_s = []
            
            for i in range(self.N):
                print(" i = ", i)
                M_omega_dot_s_i = []
                
                for j in range(3*self.N):
                    print("  j = ", j)
                    M_omega_dot_s_ij = []
                    
                    for k in range(3*self.N):
                        print("   k = ", k)
                        
                        M_omega_dot_s_ij.append(M_omega_dot(i, j, k, s))
                    M_omega_dot_s_i.append(M_omega_dot_s_ij)
                M_omega_dot_s_diff_by_s.append(sy.Matrix(M_omega_dot_s_i))
            M_omega_dot_s.append(M_omega_dot_s_diff_by_s)
        
        
        self.M_omega_dot_s = M_omega_dot_s
        
        
        
        
        
        
        print("M_omega_dot done!")



if __name__ == "__main__":
    
    
    hoge = Dynamics(3)