import numpy as np
from math import sin, cos, sqrt
import matplotlib.pyplot as plt
import matplotlib.animation as anm
from scipy.integrate import solve_ivp

from kinematics import KinematicsOfOneSection
#from dynamics import Dynamics
from controller import PD_FeedbackLinearization_Controller



class Simulator:
    """シミュレーター"""
    
    sol = None
    
    def __init__(
        self, TIME_SPAN=None, TIME_INTERVAL=None,
        pd=None, pd_dot=None, pd_dot_dot=None
    ):
        
        if TIME_SPAN is not None:
            self.TIME_SPAN = TIME_SPAN
            self.TIME_INTERVAL = TIME_INTERVAL
        
        self.Kd = 200
        self.Kp = 10000
        
        self.kinematics = KinematicsOfOneSection()
        
        
        if pd is not None and pd_dot is not None and pd_dot is not None:
            self.set_desired_position(pd, pd_dot, pd_dot_dot)
    
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
        
        self.x_range = (x_mid-max_range, x_mid+max_range)
        self.y_range = (y_mid-max_range, y_mid+max_range)
        self.z_range = (z_mid-max_range, z_mid+max_range)
    
    
    
    def set_desired_position(self, pd, pd_dot, pd_dot_dot):
        """目標位置の関数をセット"""
        
        self.pd = pd
        self.pd_dot = pd_dot
        self.pd_dot_dot = pd_dot_dot
        
        



    def calc_q_dot_dot(self, q, q_dot, p, p_dot, J, J_dot, pd, pd_dot, pd_dot_dot):
        """アクチュエータ空間上の加速度を計算"""
        #print(np.linalg.det(J))  # これが発散?
        #print("error = ", np.linalg.norm(p - pd))
        if np.linalg.det(J) < 1e-8:
            print("特異姿勢!")
            #return np.zeros((3, 1))
            z = np.linalg.pinv(J) @ (pd_dot_dot - self.Kd*(p_dot - pd_dot) - self.Kp*(p - pd) - J_dot @ q_dot)
        else:
            z = np.linalg.pinv(J) @ (pd_dot_dot - self.Kd*(p_dot - pd_dot) - self.Kp*(p - pd) - J_dot @ q_dot)
        
        return z


    def state_dot(self, t, state):
        print("t = ", t)
        q = np.array([state[:3]]).T
        q_dot = np.array([state[3:]]).T
        
        # p = self.kinematics.mapping_from_actuator_to_task_p(q, xi=1)
        # J = self.kinematics.jacobian_dpdq(q, xi=1)
        
        p = self.kinematics.linearized_mapping_from_actuator_to_task_p(q, xi=1)
        J = self.kinematics.linearized_jacobian_dpdq(q, xi=1)
        J_dot = self.kinematics.linearized_jacobian_dpdq_dot(q, q_dot, xi=1)
        
        p_dot = J @ q_dot
        
        q_dot_dot = self.calc_q_dot_dot(
            q, q_dot,
            p, p_dot, J, J_dot,
            pd = self.pd(t),
            pd_dot = self.pd_dot(t),
            pd_dot_dot = self.pd_dot_dot(t),
        )
        
        z = np.concatenate([q_dot, q_dot_dot])
        
        return np.ravel(z)



    def run_simulation(self,):
        """動力学なしで軌道追従をシミュレーション"""
        
        q_init = np.array([[-0.1, -0.1, -0.1]]).T
        dq_init = np.zeros((3, 1))
        
        state_init = np.concatenate([q_init, dq_init])
        
        self.sol = solve_ivp(
            fun = self.state_dot,
            t_span = (0, self.TIME_SPAN),
            y0 = np.ravel(state_init),
            t_eval = np.arange(0, self.TIME_SPAN, self.TIME_INTERVAL)
        )
        
        if self.sol.status == -1:
            print("失敗...")
        elif self.sol.status == 0:
            print("成功!!!")
        else:
            print("途中で終了")
    
    
    
    def plot_actuator_data(self,):
        """基本的なものをプロット"""
        
        if self.sol is None:
            return
        
        fig = plt.figure(figsize=(10, 5))
        ax = fig.add_subplot(1, 2, 1)
        ax.plot(self.sol.t, self.sol.y[0], label = "l1")
        ax.plot(self.sol.t, self.sol.y[1], label = "l2")
        ax.plot(self.sol.t, self.sol.y[2], label = "l3")
        ax.set_xlabel("time [s]")
        ax.legend()
        ax.set_xlim(0, self.sol.t[-1])
        ax.grid()
        
        ax2 = fig.add_subplot(1, 2, 2)
        ax2.plot(self.sol.t, self.sol.y[3], label = "l1_dot")
        ax2.plot(self.sol.t, self.sol.y[4], label = "l2_dot")
        ax2.plot(self.sol.t, self.sol.y[5], label = "l3_dot")
        ax2.set_xlabel("time [s]")
        ax2.legend()
        ax2.set_xlim(0, self.sol.t[-1])
        ax2.grid()
        
        fig.savefig("misc/softrobot.png")
        
        plt.show()
    
    
    def plot_all(self,):
        """全部プロット"""
        
        _q_data = np.array([self.sol.y[0], self.sol.y[1], self.sol.y[2]])
        _q_dot_data = np.array([self.sol.y[3], self.sol.y[4], self.sol.y[5]])
        
        q_data = np.split(_q_data, len(self.sol.y[0]), axis=1)
        q_dot_data = np.split(_q_dot_data, len(self.sol.t), axis=1)
        
        
        ee_p_data = np.concatenate([
            self.kinematics.linearized_mapping_from_actuator_to_task_p(q, xi=1) for q in q_data
        ], axis=1)
        ee_p_dot_data = np.concatenate([
            self.kinematics.linearized_jacobian_dpdq(q, xi=1) @ q_dot for (q, q_dot) in zip(q_data, q_dot_data)
        ], axis=1)
        
        pd_data = np.concatenate(
            [self.pd(t) for t in self.sol.t], axis=1
        )
        pd_dot_data = np.concatenate(
            [self.pd_dot(t) for t in self.sol.t], axis=1
        )
        
        error = ee_p_data - pd_data
        error_dot = ee_p_dot_data - pd_dot_data
        
        error_norm = np.linalg.norm(error, axis=0)
        error_dot_norm = np.linalg.norm(error_dot, axis=0)
        
        
        fig = plt.figure(figsize=(15, 10))
        ax = fig.add_subplot(2, 2, 1)
        ax.plot(self.sol.t, self.sol.y[0], label = "l1")
        ax.plot(self.sol.t, self.sol.y[1], label = "l2")
        ax.plot(self.sol.t, self.sol.y[2], label = "l3")
        ax.set_xlabel("time [s]")
        ax.set_ylabel("[m]")
        ax.legend()
        ax.set_xlim(0, self.sol.t[-1])
        ax.grid()
        
        ax2 = fig.add_subplot(2, 2, 2)
        ax2.plot(self.sol.t, self.sol.y[3], label = "l1_dot")
        ax2.plot(self.sol.t, self.sol.y[4], label = "l2_dot")
        ax2.plot(self.sol.t, self.sol.y[5], label = "l3_dot")
        ax2.set_xlabel("time [s]")
        ax2.set_ylabel("[m/s]")
        ax2.legend()
        ax2.set_xlim(0, self.sol.t[-1])
        ax2.grid()
        
        ax3 = fig.add_subplot(2, 2, 3)
        ax3.plot(self.sol.t, error_norm, label = "position error")
        ax3.set_xlabel("time [s]")
        ax3.set_ylabel("[m]")
        ax3.legend()
        ax3.set_xlim(0, self.sol.t[-1])
        ax3.grid()
        
        ax4 = fig.add_subplot(2, 2, 4)
        ax4.plot(self.sol.t, error_dot_norm, label = "verosity error")
        ax4.set_xlabel("time [s]")
        ax4.set_ylabel("[m/s]")
        ax4.legend()
        ax4.set_xlim(0, self.sol.t[-1])
        ax4.grid()
        
        fig.savefig("misc/softrobot_all.png")
        
        plt.show()
    
    
    
    def plot_arm(self, ax, q,):
        """axにアームをプロット"""
        
        ps = self.kinematics.calc_all_task_ps(q)
        #print(ps)
        ax.grid(True)
        ax.set_xlabel('X[m]')
        ax.set_ylabel('Y[m]')
        ax.set_zlabel('Z[m]')

        ## 三軸のスケールを揃える
        # ax.set_xlim(self.x_range)
        # ax.set_ylim(self.y_range)
        # ax.set_zlim(self.z_range)
        ax.set_box_aspect((1,1,1))
        
        ps = self.kinematics.calc_all_task_ps(q)
        
        ax.plot(ps[:, 0], ps[:, 1], ps[:, 2], label="arm", marker="o")
        ax.scatter(ps[0, 0], ps[0, 1], ps[0, 2], label="w0")
        ax.scatter(ps[-1, 0], ps[-1, 1], ps[-1, 2], label="end-effector")
        ax.legend()
        
        return
    
    
    def plot_test(self, q):
        """チェック用"""
        
        fig = plt.figure()
        ax = fig.add_subplot(projection = '3d')
        
        self.plot_arm(ax, q)
        
        
        plt.show()
        
    
    
    def make_animation(self,):
        """アニメーションで挙動確認"""
        
        # if self.sol.status == -1:
        #     return
        
        # まずはデータ作成
        xd_data = np.concatenate(
            [self.pd(t).T for t in self.sol.t]
        )
        
        
        fig = plt.figure()
        ax = fig.add_subplot(projection = '3d')

        def update(i):
            """アップデート関数"""
            ax.cla()
            
            
            self.plot_arm(
                ax = ax,
                q = np.array([[
                self.sol.y[0][i],
                self.sol.y[1][i],
                self.sol.y[2][i],
                ]]).T
            )
            
            pd = self.pd(self.sol.t[i])
            ax.scatter([pd[0,0]], [pd[1,0]], [pd[2,0]], label="temp xd", marker="*", color="r")
            ax.plot(xd_data[:, 0], xd_data[:, 1], xd_data[:, 2], label="xd line")
            
            ax.legend()
            
            
            ax.text(0, 0, 0, str(self.sol.t[i]) + "[s]")



        ani = anm.FuncAnimation(
            fig = fig,
            func = update,
            frames = len(self.sol.t),
            interval = self.TIME_INTERVAL * 0.001,
            #blit=True
        )
        
        ani.save(
            filename = "misc/softrobot.gif",
            fps = 1 / self.TIME_INTERVAL,
            writer='pillow'
        )
        
        
        plt.show()



if __name__ == "__main__":
    
    print("hoge!!")