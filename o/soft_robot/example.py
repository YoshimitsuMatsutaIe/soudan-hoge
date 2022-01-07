"""シミュレーション"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.animation as anm



# import N_is_1.Theta0 as Theta0
# import N_is_1.Phi0 as Phi0

# import N_is_1.M0_0 as M0_0
# import N_is_1.M0_1 as M0_1
# import N_is_1.M0_2 as M0_2
# import N_is_1.M1_0 as M1_0
# import N_is_1.M1_1 as M1_1
# import N_is_1.M1_2 as M1_2
# import N_is_1.M2_0 as M2_0
# import N_is_1.M2_1 as M2_1
# import N_is_1.M2_2 as M2_2

# import N_is_1.G0 as G0
# import N_is_1.G1 as G1
# import N_is_1.G2 as G2

# import N_is_1.C0_0 as C0_0
# import N_is_1.C0_1 as C0_1
# import N_is_1.C0_2 as C0_2
# import N_is_1.C1_0 as C1_0
# import N_is_1.C1_1 as C1_1
# import N_is_1.C1_2 as C1_2
# import N_is_1.C2_0 as C2_0
# import N_is_1.C2_1 as C2_1
# import N_is_1.C2_2 as C2_2


import derived.ikko_dake.eqs.numpy_style.Theta as Theta0
import derived.ikko_dake.eqs.numpy_style.Phi as Phi0

import derived.ikko_dake.eqs.numpy_style.M0_0 as M0_0
import derived.ikko_dake.eqs.numpy_style.M0_1 as M0_1
import derived.ikko_dake.eqs.numpy_style.M0_2 as M0_2
import derived.ikko_dake.eqs.numpy_style.M1_0 as M1_0
import derived.ikko_dake.eqs.numpy_style.M1_1 as M1_1
import derived.ikko_dake.eqs.numpy_style.M1_2 as M1_2
import derived.ikko_dake.eqs.numpy_style.M2_0 as M2_0
import derived.ikko_dake.eqs.numpy_style.M2_1 as M2_1
import derived.ikko_dake.eqs.numpy_style.M2_2 as M2_2

import derived.ikko_dake.eqs.numpy_style.G0 as G0
import derived.ikko_dake.eqs.numpy_style.G1 as G1
import derived.ikko_dake.eqs.numpy_style.G2 as G2

import derived.ikko_dake.eqs.numpy_style.C0_0 as C0_0
import derived.ikko_dake.eqs.numpy_style.C0_1 as C0_1
import derived.ikko_dake.eqs.numpy_style.C0_2 as C0_2
import derived.ikko_dake.eqs.numpy_style.C1_0 as C1_0
import derived.ikko_dake.eqs.numpy_style.C1_1 as C1_1
import derived.ikko_dake.eqs.numpy_style.C1_2 as C1_2
import derived.ikko_dake.eqs.numpy_style.C2_0 as C2_0
import derived.ikko_dake.eqs.numpy_style.C2_1 as C2_1
import derived.ikko_dake.eqs.numpy_style.C2_2 as C2_2



def C(q_large, xi_large, q_dot_large):
    m_0_0 = C0_0.f(q_large, xi_large, q_dot_large)
    m_0_1 = C0_1.f(q_large, xi_large, q_dot_large)
    m_0_2 = C0_2.f(q_large, xi_large, q_dot_large)
    m_1_0 = C1_0.f(q_large, xi_large, q_dot_large)
    m_1_1 = C1_1.f(q_large, xi_large, q_dot_large)
    m_1_2 = C1_2.f(q_large, xi_large, q_dot_large)
    m_2_0 = C2_0.f(q_large, xi_large, q_dot_large)
    m_2_1 = C2_1.f(q_large, xi_large, q_dot_large)
    m_2_2 = C2_2.f(q_large, xi_large, q_dot_large)
    
    return np.array([
        [m_0_0, m_0_1, m_0_2],
        [m_1_0, m_1_1, m_1_2],
        [m_2_0, m_2_1, m_2_2],
    ])




def M(q_large, xi_large, q_dot_large):
    """wow"""
    m_0_0 = M0_0.f(q_large, xi_large, q_dot_large)
    m_0_1 = M0_1.f(q_large, xi_large, q_dot_large)
    m_0_2 = M0_2.f(q_large, xi_large, q_dot_large)
    m_1_0 = M1_0.f(q_large, xi_large, q_dot_large)
    m_1_1 = M1_1.f(q_large, xi_large, q_dot_large)
    m_1_2 = M1_2.f(q_large, xi_large, q_dot_large)
    m_2_0 = M2_0.f(q_large, xi_large, q_dot_large)
    m_2_1 = M2_1.f(q_large, xi_large, q_dot_large)
    m_2_2 = M2_2.f(q_large, xi_large, q_dot_large)
    
    return np.array([
        [m_0_0, m_0_1, m_0_2],
        [m_1_0, m_1_1, m_1_2],
        [m_2_0, m_2_1, m_2_2],
    ])



def G(q_large, xi_large, q_dot_large):
    g0 = G0.f(q_large, xi_large, q_dot_large)
    g1 = G1.f(q_large, xi_large, q_dot_large)
    g2 = G2.f(q_large, xi_large, q_dot_large)
    
    return np.array([
        [g0],
        [g1],
        [g2],
    ])




def X_dot(t, X):
    """状態方程式"""
    
    X = X.reshape(6, 1)
    q = X[:3, :]
    q_dot = X[3:, :]
    
    
    tau = np.array([[0,0,0]]).T  # 入力
    
    xi = np.array([[1,1,1]]).T
    
    inv_M = np.linalg.inv(M(q, q_dot, xi))
    
    x1_dot = q_dot
    x2_dot = -inv_M @ C(q, q_dot, xi) @ q_dot +\
        inv_M @ (tau - G(q, q_dot, xi))
    #print(np.ravel(np.concatenate([x1_dot, x2_dot])))
    #print(t)
    return np.ravel(np.concatenate([x1_dot, x2_dot]))


def run_N_is_1():
    
    X = np.array([0.01, 0.02, 0, 0, 0, 0])
    
    TIME_SPAN = 10
    TIME_INTERVAL = 0.01
    
    
    sol = solve_ivp(
        fun=X_dot,
        t_span=(0, TIME_SPAN),
        t_eval=np.arange(0, TIME_SPAN, TIME_INTERVAL),
        y0 = X,
    )
    
    

    fig = plt.figure(figsize=(12, 8))
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
        xi = np.array([[1]])
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
    ax.set_ylabel("[m]")
    ax.legend()
    ax.grid()
    
    ax = fig.add_subplot(2, 2, 4)
    ax.plot(sol.t, ee_x_dot, label="x_dot")
    ax.plot(sol.t, ee_y_dot, label="y_dot")
    ax.plot(sol.t, ee_z_dot, label="z_dot")
    ax.set_xlim(0, max(sol.t))
    ax.set_xlabel("time [s]")
    ax.set_ylabel("[m/s]")
    ax.legend()
    ax.grid()
    
    
    plt.savefig("temp.png")
    


    # アニメ
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')

    x_max = 0.1
    x_min = -0.1
    y_max = 0.1
    y_min = -0.1
    z_max = 0.2
    z_min = 0
    max_range = max(x_max-x_min, y_max-y_min, z_max-z_min)*0.5
    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2
    z_mid = (z_max + z_min) / 2


    def update(i):
        """アップデート関数"""
        ax.cla()
        
        ax.grid(True)
        ax.set_xlabel('X[m]')
        ax.set_ylabel('Y[m]')
        ax.set_zlabel('Z[m]')

        ## 三軸のスケールを揃える

        ax.set_xlim(x_mid-max_range, x_mid+max_range)
        ax.set_ylim(y_mid-max_range, y_mid+max_range)
        ax.set_zlim(z_mid-max_range, z_mid+max_range)
        ax.set_box_aspect((1,1,1))
        
        
        q = np.array([[
            sol.y[0][i],
            sol.y[1][i],
            sol.y[2][i],
        ]]).T
        
        q_dot = np.array([[
            sol.y[3][i],
            sol.y[4][i],
            sol.y[5][i],
        ]]).T
        
        
        Xs = np.concatenate(
            [Phi0.f(q, np.array([[j]]), q_dot).T for j in np.linspace(0, 1, 10)]
        )
        ax.plot(Xs[:, 0], Xs[:, 1], Xs[:, 2], label="arm", marker="o")
        
        ax.scatter(ee_x[i], ee_y[i], ee_z[i], label="ee")
        
        #xd = self.xd(self.sol.t[i])
        
        
        #ax.scatter([xd[0,0]], [xd[1,0]], [xd[2,0]], label="temp xd")
        #ax.plot(xd_data[:, 0], xd_data[:, 1], xd_data[:, 2], label="xd line")
        
        ax.legend()
        
        
        ax.text(0, 0, 0, str(TIME_INTERVAL * i) + "[s]")



    ani = anm.FuncAnimation(
        fig = fig, 
        func = update, 
        frames = int(TIME_SPAN / TIME_INTERVAL)-1,
        interval = TIME_INTERVAL * 0.001
    )
    
    ani.save(
        filename = "softrobot.gif", 
        fps = 1 / TIME_INTERVAL, 
        writer='pillow'
    )
    
    
    plt.show()


def tes():
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')

    x_max = 0.1
    x_min = -0.1
    y_max = 0.1
    y_min = -0.1
    z_max = 0.2
    z_min = 0
    max_range = max(x_max-x_min, y_max-y_min, z_max-z_min)*0.5
    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2
    z_mid = (z_max + z_min) / 2



    ax.grid(True)
    ax.set_xlabel('X[m]')
    ax.set_ylabel('Y[m]')
    ax.set_zlabel('Z[m]')

    ## 三軸のスケールを揃える

    # ax.set_xlim(x_mid-max_range, x_mid+max_range)
    # ax.set_ylim(y_mid-max_range, y_mid+max_range)
    # ax.set_zlim(z_mid-max_range, z_mid+max_range)
    # ax.set_box_aspect((1,1,1))
    
    
    q = np.array([[
        0,
        0,
        0,
    ]]).T
    
    q_dot = q
    
    
    Xs = np.concatenate(
        [Phi0.f(q, np.array([[i]]), q_dot).T for i in np.linspace(0, 1, 10)]
    )

    ax.plot(Xs[:, 0], Xs[:, 1], Xs[:, 2], label="arm", marker="o")
    

    
    #xd = self.xd(self.sol.t[i])
    ax.legend()
    plt.savefig("tes.png")


if __name__ == "__main__":
    
    #tes()

    run_N_is_1()
