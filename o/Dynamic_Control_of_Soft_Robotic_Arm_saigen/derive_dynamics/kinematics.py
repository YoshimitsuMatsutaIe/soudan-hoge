import sympy as sy
from sympy import sqrt


class Local:
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
    
    
    def As(self, q):
        l1 = q[0,0]
        l2 = q[1,0]
        l3 = q[2,0]
        A1 = l1**2 + l2**2 + l3**2 - l1*l2 - l1*l3 - l2*l3
        A2 = 2*l1 - l2 - l3
        A3 = l2 - l3
        A4 = 3*self.L0 + l1 + l2 + l3
        return A1, A2, A3, A4
    
    def _P1(self, A1, A2, A3, A4, xi):
        return -(A2 * A1**4 * A4 * xi**10) / ((self.c1 * self.r**9)) + \
            (A2 * A1**3 * A4 * xi**8) / (self.c2 * self.r**7) - \
                (A2 * A1**2 * A4 * xi**6) / (self.c3 * self.r**5) + \
                    (A2 * A1 * A4 * xi**4) / (self.c4 * self.r**3) - \
                        (A2 * A4 * xi**2) / (self.c5 * self.r)
    
    
    def _P2(self, A1, A2, A3, A4, xi):
        return -(self.sq3 * A4 * A3 * A1**4 * xi**10) / (self.c1 * self.r**9) + \
            (self.sq3 * A4 * A3 * A1**3 * xi**8) / (self.c2 * self.r**7) - \
                (self.sq3 * A4 * A3 * A1**2 * xi**6) / (self.c3 * self.r**5) + \
                    (self.sq3 * A4 * A1 * A2 * xi**4) / (self.c4 * self.r**3) - \
                        (self.sq3 * A4 * A3 * xi**2) / (self.c5 * self.r)
    
    
    def _P3(self, A1, A2, A3, A4, xi):
        return (2 * A1**4 * A4 * xi**9) / (self.c6 * self.r**8) - \
            (4 * A1**3 * A4 * xi**7) / (self.c7 * self.r**6) + \
                (2 * A1**2 * A4 * xi**5) / (self.c8 * self.r**4) - \
                    (2 * A1 *A4 * xi**3) / (self.c9 * self.r**2) + \
                        (A4 * xi) / 3
    
    
    
    
    def P(self, q, xi):
        """線形化されたアクチュエータ空間からタスク空間への写像
        
        順運動学
        """

        A1, A2, A3, A4 = self.As(q)

        return sy.Matrix([[
            self._P1(A1, A2, A3, A4, xi),
            self._P2(A1, A2, A3, A4, xi),
            self._P3(A1, A2, A3, A4, xi)
        ]]).T
    
    
    def _R11(self, A1, A2, A3, A4, xi):
        return 1 - (A2**2 * A1**4 * xi**10) / (self.c1 * self.r**10) + \
            (A2**2 * A1**3 * xi**8) / (self.c1 * self.r**8) - \
                (A2**2 * A1**2 * xi**6) / (self.c3 * self.r**6) + \
                    (A1 * A2**2 * xi**4) / (self.c4 * self.r**4) - \
                        (A2**2 * xi**2) / (self.c5 * self.r**2)

    def _R12(self, A1, A2, A3, A4, xi):
        return (self.sq3 * A2 * A3 * A1**4 * xi**10) / (self.c1 * self.r**10) + \
            (self.sq3 * A2 * A3 * A1**3 * xi**8) / (self.c2 * self.r**8) - \
                (self.sq3 * A2 * A3 * A1**2 * xi**6) / (self.c3 * self.r**6) + \
                    (self.sq3 * A2 * A3 * A1 * xi**4) / (self.c4 * self.r**4) - \
                        (self.sq3 * A2 * A3 * xi**2) / (self.c5 * self.r**2)

    def _R13(self, A1, A2, A3, A4, xi):
        return -(2 * A2 * A1**4 * xi**9) / (self.c6 * self.r**9) + \
            (4 * A2 * A1**3 * xi**7) / (self.c7 * self.r**7) - \
                (2 * A2 * A1**2 * xi**5) / (self.c8 * self.r**5) + \
                    (2 * A2 * A1 * xi**3) / (self.c9 * self.r**3) - \
                        (A2 * xi) / (3 * self.r)

    def _R21(self, A1, A2, A3, A4, xi):
        return self._R12(A1, A2, A3, A4, xi)

    def _R22(self, A1, A2, A3, A4, xi):
        return 1 - (A3**2 * A1**4 * xi**10) / (self.c10 * self.r**10) + \
            (A3**2 * A1**3 * xi**8) / (self.c11 * self.r**8) - \
                (A3**2 * A1**2 * xi**6) / (self.c12 * self.r**6) + \
                    (A3**2 * A1 * xi**4) / (self.c13 * self.r**4) - \
                        (A3**2 * xi**2) / (6 * self.r**2)

    def _R23(self, A1, A2, A3, A4, xi):
        return -(2*self.sq3 * A3 * A1**4 * xi**9) / (self.c6 * self.r**9) + \
            (4*self.sq3 * A3 * A1**3 * xi**7) / (self.c7 * self.r**7) - \
                (2*self.sq3 * A3 * A1**2 * xi**5) / (self.c8 * self.r**5) + \
                    (2*self.sq3 * A3 * A1 * xi**3) / (self.c9 * self.r**3) - \
                        (self.sq3 * A3 * xi) / (3 * self.r)

    def _R31(self, A1, A2, A3, A4, xi):
        return -self._R13(A1, A2, A3, A4, xi)

    def _R32(self, A1, A2, A3, A4, xi):
        return -self._R23(A1, A2, A3, A4, xi)

    def _R33(self, A1, A2, A3, A4, xi):
        return 1 - (2 * xi**2 * A1) / (9 * self.r**2) + \
            (2 * xi**4 * A1**2) / (self.c14 * self.r**4) - \
                (4 * xi**6 * A1**3) / (self.c3 * self.r**6) + \
                    (2 * xi**8 * A1**4) / (self.c15 * self.r**8) - \
                        (4 * xi**10 * A1**5) / (self.c1 * self.r**10)

    def R(self, q, xi):
        """線形化された回転行列"""
        
        A1, A2, A3, A4 = self.As(q)
        
        return sy.Matrix([
            [self._R11(A1, A2, A3, A4, xi), self._R12(A1, A2, A3, A4, xi), self._R13(A1, A2, A3, A4, xi)],
            [self._R21(A1, A2, A3, A4, xi), self._R22(A1, A2, A3, A4, xi), self._R23(A1, A2, A3, A4, xi)],
            [self._R31(A1, A2, A3, A4, xi), self._R32(A1, A2, A3, A4, xi), self._R33(A1, A2, A3, A4, xi)],
        ])


    def MHTM(self, q, xi):
        """モーダル同時変換行列
        
        線形化されたHomogeneous Transformation Matrix
        """
        A1, A2, A3, A4 = self.As(q)
        
        return sy.Matrix([
            [self._R11(A1, A2, A3, A4, xi), self._R12(A1, A2, A3, A4, xi), self._R13(A1, A2, A3, A4, xi), self._P1(A1, A2, A3, A4, xi)],
            [self._R21(A1, A2, A3, A4, xi), self._R22(A1, A2, A3, A4, xi), self._R23(A1, A2, A3, A4, xi), self._P2(A1, A2, A3, A4, xi)],
            [self._R31(A1, A2, A3, A4, xi), self._R32(A1, A2, A3, A4, xi), self._R33(A1, A2, A3, A4, xi), self._P3(A1, A2, A3, A4, xi)],
            [0, 0, 0, 1]
        ])



class Global(Local):
    
    
    
    def __init__(self, N):
        
        self.N = N
        # self.q_large = q_large
        # self.xi_large = xi_large
        
        self.q_large = sy.Matrix(sy.MatrixSymbol('q_large', 3*N, 1))
        self.xi_large = sy.Matrix(sy.MatrixSymbol('xi_large', N, 1))
        self.set_local()
        self.set_global()
        self.set_J_OMEGA()
        self.set_J_v()
    
    
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
        """J_OMEGAを計算"""
        
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
        
        def J_v_ij(i, j):
            if j < i:
                return self.R_s[i].T * \
                    (J_v_s[i-1][:, j:j+1] + (self.J_OMEGA_s[i][:, 3*j:3*j+3] * self.P_s[i]))
            elif i == j:
                return self.R_s[i].T * sy.diff(self.P_s[i], self.q_large[j, 0])
            else:
                return sy.zeros(3, 1)
        
        J_v_s = []
        for i in range(self.N):
            #print("i = ", i)
            J_v_s_i = []
            for j in range(self.N):
                #print("j = ", j)
                J_v_s_i.append(J_v_ij(i, j))
            J_v_s.append(sy.Matrix([J_v_s_i]))
        
        self.J_v_s = J_v_s
    
    
    def set_H_OMEGA(self,):
        
        
        
        
        pass
    
    
    def set_H_v(self,):
        
        
        
        pass


if __name__ == "__main__":
    

    N = 3

    
    hoge = Global(
        # sy.Matrix(q_large), 
        # sy.Matrix(xi_large),
        N
    )
    #print(hoge.J_v_s)
    
