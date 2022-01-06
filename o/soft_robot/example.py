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
    
    
    tau = np.array([[100,0,0]]).T  # 入力
    
    xi = np.array([[1,1,1]]).T
    
    inv_M = np.linalg.inv(M.f(q, q_dot, xi))
    
    x1_dot = q_dot
    x2_dot = - inv_M @ C.f(q, q_dot, xi) @ q_dot +\
        inv_M @ (tau - G.f(q, q_dot, xi))
    #print(np.ravel(np.concatenate([x1_dot, x2_dot])))
    #print(t)
    return np.ravel(np.concatenate([x1_dot, x2_dot]))


def run_N_is_1():
    
    X = np.array([0.01, 0.02, 0, 0, 0, 0])
    
    
    sol = solve_ivp(
        fun=X_dot,
        t_span=(0, 10),
        y0 = X,
    )
    
    

    fig = plt.figure()
    ax = fig.add_subplot(2, 2, 1)
    ax.plot(sol.t, sol.y[0], label="l1")
    ax.plot(sol.t, sol.y[1], label="l2")
    ax.plot(sol.t, sol.y[2], label="l3")
    ax.set_xlim(0, max(sol.t))
    ax.set_xlabel("time [s]")
    ax.set_ylabel("length [m]")
    ax.legend()
    ax.grid()
    
    ax = fig.add_subplot(2, 2, 2)
    ax.plot(sol.t, sol.y[3], label="l1_dot")
    ax.plot(sol.t, sol.y[4], label="l2_dot")
    ax.plot(sol.t, sol.y[5], label="l3_dot")
    ax.set_xlim(0, max(sol.t))
    ax.set_xlabel("time [s]")
    ax.set_ylabel("[m/s]")
    ax.legend()
    ax.grid()
    
    
    
    ee_x, ee_y, ee_z = [], [], []
    ee_x_dot, ee_y_dot, ee_z_dot = [], [], []
    for i in range(len(sol.t)):
        q = np.array([[sol.y[0][i], sol.y[1][i], sol.y[2][i]]]).T
        q_dot = np.array([[sol.y[3][i], sol.y[4][i], sol.y[5][i]]]).T
        xi = np.array([[1,1,1]]).T
        ee = Phi0.f(q, xi, q_dot)
        ee_x.append(ee[0,0])
        ee_y.append(ee[1,0])
        ee_z.append(ee[2,0])
        
        ee_dot = Theta0.f(q, xi, q_dot) @ q_dot
        ee_x_dot.append(ee_dot[0,0])
        ee_y_dot.append(ee_dot[1,0])
        ee_z_dot.append(ee_dot[2,0])
        
    
    ax = fig.add_subplot(2, 2, 3)
    ax.plot(sol.t, ee_x, label="x")
    ax.plot(sol.t, ee_y, label="y")
    ax.plot(sol.t, ee_z, label="z")
    ax.set_xlim(0, max(sol.t))
    ax.set_xlabel("time [s]")
    ax.set_ylabel("m")
    ax.legend()
    ax.grid()
    
    ax = fig.add_subplot(2, 2, 4)
    ax.plot(sol.t, ee_x_dot, label="x_dot")
    ax.plot(sol.t, ee_y_dot, label="y_dot")
    ax.plot(sol.t, ee_z_dot, label="z_dot")
    ax.set_xlim(0, max(sol.t))
    ax.set_xlabel("time [s]")
    ax.set_ylabel("m/s")
    ax.legend()
    ax.grid()
    
    
    plt.savefig("temp.png")
    
    plt.show()




if __name__ == "__main__":
    

    run_N_is_1()
