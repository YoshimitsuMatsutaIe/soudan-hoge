"""ダイナミクスを導出"""

import sympy as sy
#from sympy import sqrt
import tqdm

import kinematics
import utils

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
    k_e = 1700  # 弾性係数
    
    def __init__(self, N,):
        """
        
        Parameters
        ---
        N : int
            セクションの数
        
        N = 3くらいでも計算にかなり時間がかかる
        """
        
        super().__init__(N)
        
        # self.set_M_omega()
        
        # self.set_M_omega_dot()
        
        # self.set_M_v()
        
        # self.set_M_v_dot()
        
        # self.set_M()
        
        # self.set_M_dot()
        
        # self.set_C()
        
        # self.set_G_g()
        
        # self.set_G_e()
        
        # self.set_G()





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
                A = sy.integrate(z, (self.xi_large[i, 0], 0, 1))
                A = operate_V(A).T
                
                B = operate_V(self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T)
                B.subs(self.xi_large[i, 0], 1)
                
                return self.Ixx * A * B
            
            elif 3*i <= j <= 3*i+2 and 3*i <= k <= 3*i+2:
                Ri = self.R_s[i]
                Ri_diff_j = sy.diff(Ri, self.q_large[j, 0])
                Ri_diff_k = sy.diff(Ri, self.q_large[k, 0])
                z = Ri_diff_j.T * Ri_diff_k
                A = sy.integrate(z, (self.xi_large[i, 0], 0, 1))
                A = operate_T2(A)
                
                return self.Ixx * A
            
            else:
                return 0

        
        M_omega_s = []
        for i in tqdm.tqdm(range(self.N)):
            #print("i = ", i)
            M_omega_s_i = []
            
            for j in tqdm.tqdm(range(3*self.N), leave=False):
                #print("j = ", j)
                M_omega_s_ij = []
                
                for k in tqdm.tqdm(range(3*self.N), leave=False):
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
                    A = sy.integrate(RR_T, (self.xi_large[i, 0], 0, 1))
                    A = operate_V(A).T
                    
                    B = operate_V(self.H_OMEGA_s[i-1][3*j:3*j+3, 3*s:3*s+3].T)
                    
                    Z = self.Ixx * A * B
                    return Z.subs(self.xi_large[i, 0], 1)
                
                elif 3*i <= s <= 3*i+2:
                    RR_T_diff_s = sy.diff(RR_T, self.q_large[s, 0])
                    A = sy.integrate(RR_T_diff_s, (self.xi_large[i, 0], 0, 1))
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
                    A = sy.integrate(z, (self.xi_large[i, 0], 0, 1))
                    A = operate_T2(A)
                    A = sy.diff(A, self.q_large[s, 0])
                    return self.Ixx * A
                
                else:
                    return 0
            
            else:
                return 0

        
        M_omega_dot_s = []
        
        for s in tqdm.tqdm(range(3*self.N)):
            #print("s = ", s)
            M_omega_dot_s_diff_by_s = []
            
            for i in tqdm.tqdm(range(self.N), leave=False):
                #print(" i = ", i)
                M_omega_dot_s_i = []
                
                for j in tqdm.tqdm(range(3*self.N), leave=False):
                    #print("  j = ", j)
                    M_omega_dot_s_ij = []
                    
                    for k in tqdm.tqdm(range(3*self.N), leave=False):
                        #print("   k = ", k)
                        
                        M_omega_dot_s_ij.append(M_omega_dot(i, j, k, s))
                    M_omega_dot_s_i.append(M_omega_dot_s_ij)
                M_omega_dot_s_diff_by_s.append(sy.Matrix(M_omega_dot_s_i))
            M_omega_dot_s.append(M_omega_dot_s_diff_by_s)
        
        
        self.M_omega_dot_s = M_omega_dot_s
        
        
        print("M_omega_dot done!")


    def set_M_v(self,):
        """並進運動の慣性行列をセット"""
        print("computing M_v ...")
        
        def M_v(i, j, k):
            if j < 3*i and k < 3*i:
                integrated_Pi = sy.integrate(
                    self.P_s[i],
                    (self.xi_large[i, 0], 0, 1)
                )
                
                A = self.m * self.J_v_s[i-1][:, j:j+1].T *\
                    (self.J_v_s[i-1][:, k:k+1] +\
                        self.J_OMEGA_s[i-1][:, 3*k:3*k+3] * integrated_Pi)
                
                B = self.m * integrated_Pi.T *\
                    self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T *\
                        self.J_v_s[i-1][:, k:k+1]
                
                integrated_PiPi_T = sy.integrate(
                    self.P_s[i] * self.P_s[i].T,
                    (self.xi_large[i, 0], 0, 1)
                )
                
                C = self.m *\
                    operate_V(integrated_PiPi_T).T *\
                        operate_V(
                            self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T * self.J_OMEGA_s[i-1][:, 3*k:3*k+3]
                        )
                
                return (A + B + C).subs(self.xi_large[i, 0], 1)
            
            
            elif j < 3*i and 3*i <= k <= 3*i+2:
                
                Pi_diff_k = sy.diff(self.P_s[i], self.q_large[k, 0])
                
                integrated_Pi_diff_k = sy.integrate(
                    Pi_diff_k,
                    (self.xi_large[i, 0], 0, 1)
                )
                
                A = self.m * self.J_v_s[i-1][:, j:j+1].T *\
                    integrated_Pi_diff_k
                
                integrated_PiPi_diff_k_T = sy.integrate(
                    self.P_s[i] *  Pi_diff_k.T,
                    (self.xi_large[i, 0], 0, 1)
                )
                
                B = self.m *\
                    operate_V(integrated_PiPi_diff_k_T).T *\
                        operate_V(self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T)
                
                return (A + B).subs(self.xi_large[i, 0], 1)
            
            
            elif 3*i <= j <= 3*i+2 and 3*i <= k <= 3*i+2:
                Pi_diff_j = sy.diff(self.P_s[i], self.q_large[j, 0])
                Pi_diff_k = sy.diff(self.P_s[i], self.q_large[k, 0])
                
                integrated = sy.integrate(
                    Pi_diff_j.T * Pi_diff_k,
                    (self.xi_large[i, 0], 0, 1)
                )
                A = self.m * integrated

                return A.subs(self.xi_large[i, 0], 1)
            
            else:
                return 0

        
        M_v_s = []
        for i in tqdm.tqdm(range(self.N)):
            #print("i = ", i)
            M_v_s_i = []
            
            for j in tqdm.tqdm(range(3*self.N), leave=False):
                #print(" j = ", j)
                M_v_s_ij = []
                
                for k in tqdm.tqdm(range(3*self.N), leave=False):
                    #print("  k = ", k)
                    M_v_s_ij.append(M_v(i, j, k))
                M_v_s_i.append(M_v_s_ij)
            M_v_s.append(sy.Matrix(M_v_s_i))
        
        
        self.M_v_s = M_v_s
        
        print("done computing M_v!")


    def set_M_v_dot(self,):
        """並進運動慣性行列の微分をセット"""
        print("computing M_v_dot ...")
        
        
        def M_v_dot(i, j, k, s):
            if j < 3*i and k < 3*i:
                if s < 3*i:
                    integrated_Pi = sy.integrate(
                        self.P_s[i],
                        (self.xi_large[i, 0], 0, 1)
                    )
                    
                    A = self.m * self.H_v_s[i-1][3*j:3*j+3, s:s+1].T *\
                        (self.J_v_s[i-1][:, k:k+1] + \
                            self.J_OMEGA_s[i-1][:, 3*k:3*k+3] * integrated_Pi)
                    
                    B = self.m * self.J_v_s[i-1][:, j:j+1].T *\
                        (self.H_v_s[i-1][3*k:3*k+3, s:s+1] +\
                            self.H_OMEGA_s[i-1][3*k:3*k+3, 3*s:3*s+3] * integrated_Pi)
                    
                    C = self.m * integrated_Pi.T *\
                        self.H_OMEGA_s[i-1][3*j:3*j+3, 3*s:3*s+3].T *\
                            self.J_v_s[i-1][:, k:k+1]

                    D = self.m * integrated_Pi.T *\
                        self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T *\
                            self.H_v_s[i-1][3*k:3*k+3, s:s+1]
                    
                    integrated_PiPi_T = sy.integrate(
                        self.P_s[i] * self.P_s[i].T,
                        (self.xi_large[i, 0], 0, 1)
                    )

                    E = self.m *\
                        operate_V(integrated_PiPi_T).T *\
                            operate_V(
                                self.H_OMEGA_s[i-1][3*j:3*j+3, 3*s:3*s+3].T *\
                                    self.J_OMEGA_s[i-1][:, 3*k:3*k+3] +\
                                        self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T *\
                                            self.H_OMEGA_s[i-1][3*k:3*k+3, 3*s:3*s+3]
                            )
                    
                    Z = A + B + C + E
                    return Z.subs(self.xi_large[i, 0], 1)
                
                elif 3*i <= s <= 3*i+2:
                    Pi_diff_s = sy.diff(self.P_s[i], self.q_large[s, 0])
                    integrated_Pi_diff_s = sy.integrate(
                        Pi_diff_s,
                        (self.xi_large[i, 0], 0, 1)
                    )
                    
                    PiPi_T_diff_s = sy.diff(
                        self.P_s[i] * self.P_s[i].T,
                        self.q_large[s, 0]
                    )
                    integrated_PiPi_T_diff_s = sy.integrate(
                        PiPi_T_diff_s,
                        (self.xi_large[i, 0], 0, 1)
                    )
                    
                    A = self.m * self.J_v_s[i-1][:, j:j+1].T *\
                        self.J_OMEGA_s[i-1][:, 3*k:3*k+3] *\
                            integrated_Pi_diff_s
                    
                    B = self.m * integrated_Pi_diff_s.T *\
                        self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T *\
                            self.J_v_s[i-1][:, k:k+1]
                    
                    C = self.m *\
                        operate_V(integrated_PiPi_T_diff_s).T *\
                            operate_V(
                                self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T *\
                                    self.J_OMEGA_s[i-1][:, 3*k:3*k+3]
                            )
                    
                    Z = A * B + C
                    return Z.subs(self.xi_large[i, 0], 1)
                
                else:
                    return 0
            
            
            elif j < 3*i and 3*i <= k <= 3*i+2:
                Pi_diff_k = sy.diff(self.P_s[i], self.q_large[k, 0])

                if s < 3*i:
                    
                    integrated_Pi_diff_k = sy.integrate(
                        Pi_diff_k,
                        (self.xi_large[i, 0], 0, 1)
                    )
                    
                    A = self.m * self.H_v_s[i-1][3*j:3*j+3, s:s+1].T *\
                        integrated_Pi_diff_k
                    
                    integrated_PiPi_diff_k_T = sy.integrate(
                        self.P_s[i] * Pi_diff_k.T,
                        (self.xi_large[i, 0], 0, 1)
                    )
                    
                    B = self.m *\
                        operate_V(integrated_PiPi_diff_k_T).T *\
                            operate_V(self.H_OMEGA_s[i-1][3*j:3*j+3, 3*s:3*s+3])
                    
                    Z = A + B
                    return Z.subs(self.xi_large[i, 0], 1)
                
                elif 3*i <= s <= 3*i+2:
                    Pi_diff_ks = sy.diff(Pi_diff_k, self.q_large[s, 0])
                    integrated_Pi_diff_ks = sy.integrate(
                        Pi_diff_ks,
                        (self.xi_large[i, 0], 0, 1)
                    )
                    
                    A = self.m * self.J_v_s[i-1][:, j:j+1].T *\
                        integrated_Pi_diff_ks
                    
                    PiPi_diff_k_T_diff_s = sy.diff(
                        self.P_s[i] * Pi_diff_k.T,
                        self.q_large[s, 0]
                    )
                    integrated_PiPi_diff_k_T_diff_s = sy.integrate(
                        PiPi_diff_k_T_diff_s,
                        (self.xi_large[i, 0], 0, 1)
                    )
                    
                    B = self.m *\
                        operate_V(integrated_PiPi_diff_k_T_diff_s).T *\
                            operate_V(self.J_OMEGA_s[i-1][:, 3*j:3*j+3].T)
                    
                    Z = A + B
                    return Z.subs(self.xi_large[i, 0], 1)

                else:
                    return 0
            
            elif 3*i <= j <= 3*i+2 and 3*i <= k <= 3*i+2:
                if s < 3*i:
                    return 0
                
                elif 3*i <= s <= 3*i+2:
                    Pi_diff_j = sy.diff(self.P_s[i], self.q_large[j, 0])
                    Pi_diff_k = sy.diff(self.P_s[i], self.q_large[k, 0])
                    
                    Z = self.m *\
                        sy.diff(
                            Pi_diff_j.T * Pi_diff_k,
                            self.q_large[s, 0]
                        )
                    
                    return Z.subs(self.xi_large[i, 0], 1)

                else:
                    return 0
            
            else:
                return 0

        
        M_v_dot_s = []
        
        for s in tqdm.tqdm(range(3*self.N)):
            #print("s = ", s)
            M_v_dot_s_diff_by_s = []
            
            for i in tqdm.tqdm(range(self.N), leave=False):
                #print(" i = ", i)
                M_v_dot_s_i = []
                
                for j in tqdm.tqdm(range(3*self.N), leave=False):
                    #print("  j = ", j)
                    M_v_dot_s_ij = []
                    
                    for k in tqdm.tqdm(range(3*self.N), leave=False):
                        #print("   k = ", k)
                        
                        M_v_dot_s_ij.append(M_v_dot(i, j, k, s))
                    M_v_dot_s_i.append(M_v_dot_s_ij)
                M_v_dot_s_diff_by_s.append(sy.Matrix(M_v_dot_s_i))
            M_v_dot_s.append(M_v_dot_s_diff_by_s)
        
        
        self.M_v_dot_s = M_v_dot_s
        
        
        print("done computing M_v_dot!")


    def set_M(self,):
        """慣性行列をセット"""
        print("computing M...")
        
        M = sy.zeros(3*self.N, 3*self.N)
        for i in tqdm.tqdm(range(self.N)):
            M = M + self.M_v_s[i] + self.M_omega_s[i]
        
        self.M = M

        print("done computing N!")


    def set_M_dot(self,):
        """慣性行列の各微分値をセット"""
        print("computing M_dot...")
        
        M_dot_s = []
        for s in tqdm.tqdm(range(3*self.N)):
            M_dot_s_s = sy.zeros(3*self.N, 3*self.N)
            for i in tqdm.tqdm(range(self.N), leave=False):
                M_dot_s_s += self.M_omega_dot_s[s][i] + self.M_v_dot_s[s][i]
            M_dot_s.append(M_dot_s_s)
        
        self.M_dot_s = M_dot_s
        print("done computing M_dot!")


    def christoffel_M(self, k, j, h):
        """Mのクリストッフェル記号???"""
        
        return 1/2 *\
            (self.M_dot_s[h][k, j] +\
                self.M_dot_s[j][k, h] -\
                    self.M_dot_s[k][h, j])


    def set_C(self,):
        """コリオリ・遠心力をセット"""
        print("computing C...")
        
        
        def C_ikj(i, k, j):
            w = 3 * (i+1)
            
            c = 0
            for h in range(w):
                c += self.christoffel_M(k, j, h) * self.q_dot_large[h, 0]
            
            return c
        
        
        C = sy.zeros(3*self.N, 3*self.N)
        
        for i in tqdm.tqdm(range(self.N)):
            for k in tqdm.tqdm(range(3*self.N), leave=False):
                for j in tqdm.tqdm(range(3*self.N), leave=False):
                    C[k, j] += C_ikj(i, k, j)

        self.C = C
        print("done computing C !")


    def set_G_g(self,):
        """一般化重力ベクトルをセット"""
        print("computing G_g ...")
        
        def G_g(i, j):
            if j < 3*i:
                integrated_Pi = sy.integrate(
                    self.P_s[i],
                    (self.xi_large[i, 0], 0, 1)
                )
                Z = self.m *\
                    (self.J_v_s[i-1][:, j:j+1] +\
                        self.J_OMEGA_s[i-1][:, 3*j:3*j+3] * integrated_Pi).T *\
                            self.Theta_s[i-1].T *\
                                self.g
                
                return Z.subs(self.xi_large[i, 0], 1)
            
            elif 3*i <= j <= 3*i*2:
                Pi_T_diff_j = sy.diff(
                    self.P_s[i].T,
                    self.q_large[j, 0]
                )
                integrated = sy.integrate(
                    Pi_T_diff_j,
                    (self.xi_large[i, 0], 0, 1)
                )
                
                Z = self.m * integrated * self.Theta_s[i-1].T * self.g
                return Z.subs(self.xi_large[i, 0], 1)
            
            else:
                return 0
        
        
        G_g_s = []
        for i in tqdm.tqdm(range(self.N)):
            
            G_g_s_i = []
            for j in tqdm.tqdm(range(3*self.N), leave=False):
                G_g_s_i.append([G_g(i, j)])
            G_g_s.append(sy.Matrix(G_g_s_i))
        
        self.G_g_s = G_g_s
        print("done computing G_g !")


    def set_G_e(self,):
        """弾性ポテンシャルエネルギーに関する力をセット"""
        print("computing G_e...")
        self.G_e = sy.diag([self.k_e for _ in range(3*self.N)], unpack=True) * self.q_large
        print("done computing G_e!")


    def set_G(self,):
        """G行列をセット"""
        print("computing G...")
        
        G = sy.zeros(3*self.N, 1)
        
        for i in tqdm.tqdm(range(self.N)):
            G += self.G_g_s[i]
        
        G += self.G_e
        
        self.G = G
        print("done computing G!")


if __name__ == "__main__":
    
    
    hoge = Dynamics(1)