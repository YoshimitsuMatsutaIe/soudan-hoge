import sympy as sy
from sympy import sqrt


class Base:
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
    
    
    
    
    def local_P(self, q, xi):
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

    def local_R(self, q, xi):
        """線形化された回転行列"""
        
        A1, A2, A3, A4 = self.As(q)
        
        return sy.Matrix([
            [self._R11(A1, A2, A3, A4, xi), self._R12(A1, A2, A3, A4, xi), self._R13(A1, A2, A3, A4, xi)],
            [self._R21(A1, A2, A3, A4, xi), self._R22(A1, A2, A3, A4, xi), self._R23(A1, A2, A3, A4, xi)],
            [self._R31(A1, A2, A3, A4, xi), self._R32(A1, A2, A3, A4, xi), self._R33(A1, A2, A3, A4, xi)],
        ])


    def local_MHTM(self, q, xi):
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




if __name__ == "__main__":
    
    l1, l2, l3 = sy.symbols("l1, l2, l3")
    xi = sy.Symbol("xi")
    
    q = sy.Matrix([[l1, l2, l3]]).T
    
    hoge = Base()
    print(hoge.local_MHTM(q, xi))