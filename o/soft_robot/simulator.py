"""シミュレーション"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


import N_is_1.Theta0 as Theta0
import N_is_1.Phi0 as Phi0
import N_is_1.M as M
import N_is_1.C as C
import N_is_1.G as G


def X_dot(t, X):
    """状態方程式"""
    
    X = X.reshape(6, 1)
    q = X[:3, :]
    q_dot = X[3:, :]
    
    
    tau = np.array([[0,0,0]]).T  # 入力
    
    xi = np.array([[1,1,1]]).T
    
    inv_M = np.linalg.inv(M.f(q, q_dot, xi))
    
    x1_dot = q_dot
    x2_dot = - inv_M @ C.f(q, q_dot, xi) @ q_dot +\
        inv_M @ (tau - G.f(q, q_dot, xi))
    #print(np.ravel(np.concatenate([x1_dot, x2_dot])))
    print(t)
    return np.ravel(np.concatenate([x1_dot, x2_dot]))





if __name__ == "__main__":
    
    X = np.array([0.01, 0.02, 0, 0, 0, 0])
    
    
    sol = solve_ivp(
        fun=X_dot,
        t_span=(0, 0.1),
        y0 = X,
    )
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(sol.t, sol.y[0])
    ax.plot(sol.t, sol.y[1])
    ax.plot(sol.t, sol.y[2])
    
    
    plt.show()