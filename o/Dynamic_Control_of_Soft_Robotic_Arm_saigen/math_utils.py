"""関数"""

import numpy as np


def sknew_symmetric(a):
    """^演算
    
    ベクトル -> 歪対称行列"""
    return np.array([
        [0, -a[2,0], a[1,0]],
        [a[2,0], 0, -a[0,0]],
        [-a[1,0], a[0,0], 0],
    ])


def inv_sknew_symmetric(A):
    """逆^演算
    
    歪対称行列 -> ベクトル"""
    return np.array([[A[2,1], A[0,2], A[1,0]]]).T



def T2(M):
    """主対角線上の最初の2つの要素の和"""
    return M[0,0] + M[1,1]


