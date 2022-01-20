"""運動学関係

2次のマクローリン展開
"""

import sympy as sy
from sympy import sqrt

import maclaurin_2.P_0 as P_0
import maclaurin_2.P_1 as P_1
import maclaurin_2.P_2 as P_2

import maclaurin_2.R_0_0 as R_0_0
import maclaurin_2.R_0_1 as R_0_1
import maclaurin_2.R_0_2 as R_0_2
import maclaurin_2.R_0_0 as R_1_0
import maclaurin_2.R_0_1 as R_1_1
import maclaurin_2.R_0_2 as R_1_2
import maclaurin_2.R_0_0 as R_2_0
import maclaurin_2.R_0_1 as R_2_1
import maclaurin_2.R_0_2 as R_2_2


class Local:
    """ローカル座標系関連  
    
    全セクションで同一のパラメータであると仮定  
    """
    
    def P(self, q, xi):
        """線形化されたアクチュエータ空間からタスク空間への写像
        
        順運動学
        """



        return sy.Matrix([[
            P_0.f(q, xi),
            P_1.f(q, xi),
            P_2.f(q, xi)
        ]]).T
    

    def R(self, q, xi):
        """線形化された回転行列"""
        

        return sy.Matrix([
            [R_0_0.f(q, xi), R_0_1.f(q, xi), R_0_2.f(q, xi)],
            [R_1_0.f(q, xi), R_1_1.f(q, xi), R_1_2.f(q, xi)],
            [R_2_0.f(q, xi), R_2_1.f(q, xi), R_2_2.f(q, xi)],
        ])


    def MHTM(self, q, xi):
        """モーダル同時変換行列
        
        線形化されたHomogeneous Transformation Matrix
        """
        
        return sy.Matrix([
            [R_0_0.f(q, xi), R_0_1.f(q, xi), R_0_2.f(q, xi), P_0.f(q, xi)],
            [R_1_0.f(q, xi), R_1_1.f(q, xi), R_1_2.f(q, xi), P_1.f(q, xi)],
            [R_2_0.f(q, xi), R_2_1.f(q, xi), R_2_2.f(q, xi), P_2.f(q, xi)],
            [0, 0, 0, 1]
        ])



class Global(Local):
    """位置，回転行列，ヤコビアン，ヘッシアン等のグローバル表現  
    
    全セクションで同一のパラメータであると仮定.   
    すぐ終わる.   
    """
    
    
    def __init__(self, N):
        """
        Parameters  
        ---
        N : int
            セクションの数
        """
        
        
        print("computing kinematics...")
        
        self.N = N
        # self.q_large = q_large
        # self.xi_large = xi_large
        
        self.q_large = sy.Matrix(sy.MatrixSymbol('q_large', 3*N, 1))  # 完全な関節角度ベクトル
        self.q_dot_large = sy.Matrix(sy.MatrixSymbol('q_dot_large', 3*N, 1))  # 完全な関節角速度ベクトル
        self.xi_large = sy.Matrix(sy.MatrixSymbol('xi_large', N, 1))  # 完全なスカラξベクトル
        
        self.set_local()
        self.set_global()
        self.set_J_OMEGA()
        self.set_J_v()
        self.set_H_OMEGA()
        self.set_H_v()
        
        print("done computing kinematics!")
    
    
    def set_local(self,):
        
        self.P_s = []
        self.R_s = []
        for i in range(self.N):
            q = self.q_large[i:i+3, :]
            xi = self.xi_large[i, 0]
            self.P_s.append(self.P(q, xi))
            self.R_s.append(self.R(q, xi))
    
    
    def set_global(self,):
        
        self.Theta_s = []
        self.Phi_s = []
        for i in range(self.N):
            if i == 0:
                self.Theta_s.append(self.R_s[0])
                self.Phi_s.append(self.P_s[0])
            else:
                Theta = self.Theta_s[i-1] * self.R_s[i]
                self.Theta_s.append(Theta.subs(self.xi_large[i-1, 0], 1))
                
                Phi = self.Phi_s[i-1] + self.Theta_s[i-1] * self.P_s[i]
                self.Phi_s.append(Phi.subs(self.xi_large[i-1, 0], 1))



    def set_J_OMEGA(self,):
        """角速度ヤコビアンのハット変換を計算"""
        
        def J_OMEGA_ij(i, j):
            if j <= 3*i-1:
                return self.R_s[i].T * J_OMEGA_s[i-1][:, 3*j:3*j+3] * self.R_s[i]
            elif 3*i <= j <= 3*i+2:
                return self.R_s[i].T * sy.diff(self.R_s[i], self.q_large[j, 0])
            else:
                return sy.zeros(3, 3)
        
        J_OMEGA_s = []
        for i in range(self.N):
            #print("i = ", i)
            J_OMEGAs_i = []
            for j in range(3*self.N):
                #print("j = ", j)
                J_OMEGAs_i.append(J_OMEGA_ij(i, j))
            J_OMEGA_s.append(sy.Matrix([J_OMEGAs_i]))
        
        self.J_OMEGA_s = J_OMEGA_s


    def set_J_v(self,):
        """線速度ヤコビアンをセット"""
        
        def J_v_ij(i, j):
            if j < 3*i:
                return self.R_s[i].T * \
                    (J_v_s[i-1][:, j:j+1] + (self.J_OMEGA_s[i][:, 3*j:3*j+3] * self.P_s[i]))
            elif 3*i <= j <= 3*i+2:
                return self.R_s[i].T * sy.diff(self.P_s[i], self.q_large[j, 0])
            else:
                return sy.zeros(3, 1)
        
        J_v_s = []
        for i in range(self.N):
            #print("i = ", i)
            J_v_s_i = []
            for j in range(3*self.N):
                #print("j = ", j)
                J_v_s_i.append(J_v_ij(i, j))
            J_v_s.append(sy.Matrix([J_v_s_i]))
        
        self.J_v_s = J_v_s
    
    
    def set_H_OMEGA(self,):
        """角速度ヘッシアンをセット
        
        ※本当はテンソル?  
        """
        
        
        def H_OMEGA_ijk(i, j, k):
            if j < 3*i and k < 3*i:
                H_OMEGA_prev = H_OMEGA_s[i-1][3*j:3*j+3, 3*k:3*k+3]  # テンソルから1枚剥がして持ってくる
                #print(H_OMEGA_prev.shape)
                return self.R_s[i].T * H_OMEGA_prev * self.R_s[i]

            elif j < 3*i and 3*i <= k <= 3*i+2:
                R_i_diff_k = sy.diff(self.R_s[i], self.q_large[k, 0])
                return R_i_diff_k.T * self.J_OMEGA_s[i-1][:, 3*j:3*j+3] * self.R_s[i] +\
                    self.R_s[i].T * self.J_OMEGA_s[i-1][:, 3*j:3*j+3] * R_i_diff_k
            
            elif 3*i <= j <= 3*i+2 and k < 3*i:
                return sy.zeros(3, 3)
            
            elif 3*i <= j <= 3*i+2 and 3*i <= k <= 3*i+2:
                R_i_diff_k = sy.diff(self.R_s[i], self.q_large[k, 0])
                R_i_diff_j = sy.diff(self.R_s[i], self.q_large[j, 0])
                R_i_diff_j_diff_k = sy.diff(R_i_diff_j, self.q_large[k, 0])
                return R_i_diff_k.T * R_i_diff_j + self.R_s[i].T * R_i_diff_j_diff_k
            
            else:
                return sy.zeros(3, 3)
        
        
        H_OMEGA_s = []
        for i in range(self.N):
            #print("i = ", i)
            H_OMEGA_s_i = []
            
            for j in range(3*self.N):
                #print("j = ", j)
                H_OMEGA_s_ij = []
                
                for k in range(3*self.N):
                    #print("k = ", k)
                    H_OMEGA_s_ij.append(H_OMEGA_ijk(i, j, k))
                H_OMEGA_s_i.append([sy.Matrix([H_OMEGA_s_ij])])
            
            H_OMEGA_s.append(sy.Matrix(H_OMEGA_s_i))
            #print(H_OMEGA_s[-1].shape)
        
        
        self.H_OMEGA_s = H_OMEGA_s
    
    
    def set_H_v(self,):
        """線速度ヘッシアンをセット
        
        ホントはテンソル?  
        """
        
        def H_v_ijk(i, j, k):
            if j < 3*i and k < 3*i:
                return self.R_s[i].T * \
                    (H_v_s[i-1][3*j:3*j+3, k:k+1] + self.H_OMEGA_s[i-1][3*j:3*j+3, 3*k:3*k+3] * self.P_s[i])
            
            elif j < 3*i and 3*i <= k <= 3*i+2:
                R_i_diff_k = sy.diff(self.R_s[i], self.q_large[k, 0])
                P_i_diff_k = sy.diff(self.P_s[i], self.q_large[k, 0])
                return R_i_diff_k.T *\
                    (self.J_v_s[i-1][:, j:j+1] + self.J_OMEGA_s[i-1][:, 3*j:3*j+3] * self.P_s[i]) +\
                        self.R_s[i].T * self.J_OMEGA_s[i][:, 3*j:3*j+3] * P_i_diff_k
            
            elif 3*i <= j <= 3*i+2 and k < 3*i:
                return sy.zeros(3, 1)
            
            elif 3*i <= j <= 3*i+2 and 3*i <= k <= 3*i+2:
                R_i_diff_k = sy.diff(self.R_s[i], self.q_large[k, 0])
                P_i_diff_j = sy.diff(self.P_s[i], self.q_large[j, 0])
                P_i_diff_j_diff_k = sy.diff(P_i_diff_j, self.q_large[k, 0])
                return R_i_diff_k.T * P_i_diff_j + self.R_s[i].T * P_i_diff_j_diff_k
            
            else:
                return sy.zeros(3, 1)
        
        
        H_v_s = []
        for i in range(self.N):
            #print("i = ", i)
            H_v_s_i = []
            
            for j in range(3*self.N):
                #print("j = ", j)
                H_v_s_ij = []
                
                for k in range(3*self.N):
                    #print("k = ", k)
                    H_v_s_ij.append(H_v_ijk(i, j, k))
                H_v_s_i.append([sy.Matrix([H_v_s_ij])])
            
            H_v_s.append(sy.Matrix(H_v_s_i))
            #print(H_OMEGA_s[-1].shape)
        
        
        self.H_v_s = H_v_s




if __name__ == "__main__":
    

    N = 3

    
    hoge = Global(
        # sy.Matrix(q_large), 
        # sy.Matrix(xi_large),
        N
    )
    #print(hoge.J_v_s)
    