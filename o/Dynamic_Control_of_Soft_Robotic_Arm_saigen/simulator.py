import numpy as np
from math import sin, cos, sqrt
import matplotlib.pyplot as plt
import matplotlib.animation as anm
import scipy as sp
from scipy.integrate import solve_ivp

from kinematics import Kinematics
from dynamics import Dynamics

class Simulator:
    
    sol = None
    xi_all = np.arange(0, 1, 0.01)  # xi一覧
    
    def __init__(self, TIME_SPAN, TIME_INTERVAL):
        
        self.TIME_SPAN = TIME_SPAN
        self.TIME_INTERVAL = TIME_INTERVAL
        
        self.Kd = 200
        self.Kp = 1e4
        
        self.kinematics = Kinematics()
    
    
    def xd(self, t):
        """タスク空間上の所望の位置"""
        return np.array([
            [0.1 * sin(3*t)],
            [0.1 * cos(3*t)],
            [0.147],
        ])


    def xd_dot(self, t):
        """タスク空間上の所望の測度"""
        return np.array([
            [0.1 * 3 * cos(3*t)],
            [0.1 * 3 * -sin(3*t)],
            [0],
        ])


    def xd_dot_dot(self, t):
        """タスク空間上の所望の加速度"""
        return np.array([
            [0.1 * 9 * -sin(3*t)],
            [0.1 * 9 * -cos(3*t)],
            [0],
        ])


    def calc_q_dot_dot(self, x, x_dot, J, xd, xd_dot, xd_dot_dot):
        """アクチュエータ空間上の加速度を計算"""
        #print(np.linalg.pinv(J))  # これが発散
        z = np.linalg.pinv(J) @ \
            (xd_dot_dot - self.Kd*(x_dot - xd_dot) - self.Kp*(x - xd))
        return z


    def state_dot(self, t, state):
        
        q = np.array([state[:3]]).T
        q_dot = np.array([state[3:]]).T
        
        x = self.kinematics.calc_X(q, xi=1)
        J = self.kinematics.calc_Jacobian(q, xi=1)
        x_dot = J @ q_dot
        
        q_dot_dot = self.calc_q_dot_dot(
            x, x_dot, J,
            xd = self.xd(t),
            xd_dot = self.xd_dot(t),
            xd_dot_dot = self.xd_dot_dot(t),
        )
        
        z = np.concatenate([q_dot, q_dot_dot])
        
        return np.ravel(z)



    def run_simulation(self,):
        """動力学なしで軌道追従をシミュレーション"""
        
        q_init = np.array([[0, 0, 0]]).T
        dq_init = np.zeros((3, 1))
        
        state_init = np.concatenate([q_init, dq_init])
        
        self.sol = solve_ivp(
            fun = self.state_dot,
            t_span = (0, self.TIME_SPAN),
            y0 = np.ravel(state_init),
            t_eval = np.arange(0, self.TIME_SPAN, self.TIME_INTERVAL)
        )
        print("成功したか否か", self.sol.status)

        return
    
    
    def plot_actuator_data(self,):
        """基本的なものをプロット"""
        
        if self.sol is None:
            return
        
        else:
            fig = plt.figure()
            ax = fig.add_subplot(1, 2, 1)
            ax.plot(self.sol.t, self.sol.y[0], label = "l1")
            ax.plot(self.sol.t, self.sol.y[1], label = "l2")
            ax.plot(self.sol.t, self.sol.y[2], label = "l3")
            ax.set_xlabel("time [s]")
            ax.legend()
            ax.grid()
            
            ax2 = fig.add_subplot(1, 2, 2)
            ax2.plot(self.sol.t, self.sol.y[3], label = "l1_dot")
            ax2.plot(self.sol.t, self.sol.y[4], label = "l2_dot")
            ax2.plot(self.sol.t, self.sol.y[5], label = "l3_dot")
            ax2.set_xlabel("time [s]")
            ax2.legend()
            ax2.grid()
            
            plt.show()
    
    
    def make_animation(self,):
        """アニメーションで挙動確認"""
        
        if self.sol.status == -1:
            return
        
        # まずはデータ作成
        xd_data = np.concatenate(
            [self.xd(t).T for t in self.sol.t]
        )
        
        
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
                self.sol.y[0][i],
                self.sol.y[1][i],
                self.sol.y[2][i],
            ]]).T
            
            Xs = np.concatenate(
                [self.kinematics.calc_X(q, xi).T for xi in self.xi_all]
            )
            
            
            xd = self.xd(self.sol.t[i])
            
            ax.plot(Xs[:, 0], Xs[:, 1], Xs[:, 2], label="arm")
            ax.scatter([xd[0,0]], [xd[1,0]], [xd[2,0]], label="temp xd")
            ax.plot(xd_data[:, 0], xd_data[:, 1], xd_data[:, 2], label="xd line")
            
            ax.legend()
            
            
            ax.text(0, 0, 0, str(self.TIME_INTERVAL * i) + "[s]")



        ani = anm.FuncAnimation(
            fig = fig, 
            func = update, 
            frames = int(self.TIME_SPAN / self.TIME_INTERVAL)-1,
            interval = self.TIME_INTERVAL * 0.001
        )
        
        ani.save(
            filename = "softrobot.gif", 
            fps = 1 / self.TIME_INTERVAL, 
            writer='pillow'
        )
        
        
        plt.show()



if __name__ == "__main__":
    
    hoge = Simulator(5, 0.01)
    
    hoge.run_simulation()
    hoge.plot_actuator_data()
    hoge.make_animation()
