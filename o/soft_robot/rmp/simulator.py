import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as anm
import tqdm


from kinematics import Kinematics
from differential_kinematics import DifferentialKinematics
from controller import OriginalRMPAttractor, OriginalRMPJointLimitAvoidance, pullback



class Simulator:
    """シミュレーション関係"""
    
    TIME_SPAN = 0.1  # デフォルト
    TIME_INTERVAL = 0.01  # デフォルト
    
    
    def __init__(self, N=3):
        
        self.N = N  # セクションの数
        
        self.kim = Kinematics()
        self.diff_kim = DifferentialKinematics(N)
        
        
        self.goal = np.array([[1, 1, 5]]).T
        
        
        self.q_max = np.array([[0.1] * self.N*3]).T
        self.q_min = -self.q_max
        
        self.attractor = OriginalRMPAttractor(
            max_speed = 10,
            gain = 300,
            a_damp_r = 0.05,
            sigma_W = 1,
            sigma_H = 1,
            A_damp_r = 0.01,
        )
        
        self.jlavoidance = OriginalRMPJointLimitAvoidance(
            gamma_p = 5,
            gamma_d = 1,
            lam = 1
        )
        
        

        pass
    
    
    def X_dot(self, t, X):
        """scipyで使う微分方程式"""
        print("t = ", t)
        q = X[:self.N*3].reshape(self.N*3, 1)
        q_dot = X[self.N*3:].reshape(self.N*3, 1)
        
        
        x = self.kim.Phi(self.N, q)
        J = self.diff_kim.J(self.N, q)
        x_dot = J @ q_dot
        
        root_f = np.zeros((self.N*3, 1))
        root_M = np.zeros((self.N*3, self.N*3))

        # アトラクター
        f, M = self.attractor.get_natural(x, x_dot, self.goal, np.zeros((3,1)))
        pullbacked_f, pullbacked_M = pullback(f, M, J)
        root_f += pullbacked_f
        root_M += pullbacked_M
        
        
        # ジョイント制限回避
        f, M = self.jlavoidance.get_natural(
            q = q,
            dq = q_dot,
            q_max = self.q_max,
            q_min = self.q_min
        )
        root_f += f
        root_M += M
        
        # resolve演算
        q_dot_dot = np.linalg.pinv(root_M) @ root_f
        
        return np.ravel(np.concatenate([q_dot, q_dot_dot]))
    
    
    
    def run(self, TIME_SPAN=None, TIME_INTERVAL=None):
        """シミュレーション実行"""
        
        print("シミュレーション実行中...")
        
        
        if TIME_SPAN is not None:
            self.TIME_SPAN = TIME_SPAN
        
        if TIME_INTERVAL is not None:
            self.TIME_INTERVAL = TIME_INTERVAL
        
        
        X_init = np.zeros(self.N*3*2)
        
        self.sol = integrate.solve_ivp(
            fun = self.X_dot,
            t_span = (0, self.TIME_SPAN),
            y0 = X_init,
            t_eval = np.arange(0, self.TIME_SPAN, self.TIME_INTERVAL)
        )
        
        print("シミュレーション終了!")
    
    
    def plot(self,):
        
        
        fig = plt.figure(figsize=(10, 15))
        
        # xの時系列データ作成
        e = []
        for i in range(len(self.sol.t)):
            temp_q = []
            for j in range(self.N):
                temp_q.append(self.sol.y[3*j][i])
                temp_q.append(self.sol.y[3*j+1][i])
                temp_q.append(self.sol.y[3*j+2][i])
            temp_q = np.array([temp_q]).T
            e.append(
                np.linalg.norm(self.goal - self.kim.Phi(self.N, temp_q))
            )
        
        ax = fig.add_subplot(3, 1, 1)
        ax.plot(self.sol.t, e, label="error")
        ax.set_xlabel("time [s]")
        ax.set_ylabel("[m]")
        ax.set_xlim(0, self.TIME_SPAN)
        ax.set_ylim(0, max(e))
        ax.legend()
        ax.grid()
        
        ax2 = fig.add_subplot(3, 1, 2)
        for i in range(self.N):
            ax2.plot(self.sol.t, self.sol.y[3*i], label="l" + str(i) + "_1")
            ax2.plot(self.sol.t, self.sol.y[3*i+1], label="l" + str(i) + "_2")
            ax2.plot(self.sol.t, self.sol.y[3*i+2], label="l" + str(i) + "_3")
        ax2.set_xlabel("time [s]")
        ax2.set_ylabel("length [m]")
        ax2.set_xlim(0, self.TIME_SPAN)
        ax2.legend()
        ax2.grid()
        
        ax3 = fig.add_subplot(3, 1, 3)
        for i in range(self.N):
            ax3.plot(self.sol.t, self.sol.y[3*i+3*self.N], label="dl" + str(i) + "_1")
            ax3.plot(self.sol.t, self.sol.y[3*i+1+3*self.N], label="dl" + str(i) + "_2")
            ax3.plot(self.sol.t, self.sol.y[3*i+2+3*self.N], label="dl" + str(i) + "_3")
        
        ax3.set_xlabel("time [s]")
        ax3.set_ylabel("[m/s]")
        ax3.set_xlim(0, self.TIME_SPAN)
        ax3.legend()
        ax3.grid()
        
        fig.savefig("py.png")
        plt.show()
        print("plot完了!")


    def reproduce_state(self,):
        """アクチュエータ変位から状態を再現"""
        
        print("データ復元中...")
        print("q...")
        self.q_data = []
        for i in tqdm.tqdm(range(len(self.sol.t))):
            temp = []
            for j in range(self.N):
                temp.append(self.sol.y[3*j][i])
                temp.append(self.sol.y[3*j+1][i])
                temp.append(self.sol.y[3*j+2][i])
            self.q_data.append(np.array([temp]).T)
        
        print("q_dot...")
        self.q_dot_data = []
        for i in tqdm.tqdm(range(len(self.sol.t))):
            temp = []
            for j in range(self.N):
                temp.append(self.sol.y[3*j+3*self.N][i])
                temp.append(self.sol.y[3*j+1+3*self.N][i])
                temp.append(self.sol.y[3*j+2+3*self.N][i])
            self.q_dot_data.append(np.array([temp]).T)
        
        print("x...")
        self.x_data = []
        for i in tqdm.tqdm(range(len(self.sol.t))):
            temp = []
            for j in range(self.N):
                temp.append(
                    np.concatenate(
                        [self.kim.Phi(j+1, self.q_data[i], xi) for xi in self.kim.xi_large]
                    )
                )
            self.x_data.append(temp)


        print("完了!")
        
    
    def _update(self, i):
        """アニメのフレーム作成"""
        
        t = i * self.TIME_INTERVAL
        
        
        self.ax.cla()
        
        self.ax.grid(True)
        self.ax.set_xlabel('X[m]')
        self.ax.set_ylabel('Y[m]')
        self.ax.set_zlabel('Z[m]')
        
        self.ax.set_xlim(self.xl, self.xu)
        self.ax.set_ylim(self.yl, self.yu)
        self.ax.set_zlim(self.zl, self.zu)
        self.ax.set_box_aspect((1,1,1))
        
        # 時刻表示
        self.ax.text(
            0.0, 0.0, -0.01,
            self.time_template % (i * self.TIME_INTERVAL), size = 10
        )
        
        
        
        # 目標点
        self.ax.scatter(
            self.goal[0, 0], self.goal[1, 0], self.goal[2, 0],
            s = 100, label = 'goal point', marker = '*', color = '#ff7f00', 
            alpha = 1, linewidths = 1.5, edgecolors = 'red'
        )
        
        # アーム
        for j in range(self.N):
            self.ax.plot(
                self.x_data[i][j][0,:], self.x_data[i][j][1,:], self.x_data[i][j][2,:],
                label="section" + str(j)
            )







    def make_aniation(self,):
        """アニメ作る"""
        
        print("アニメ作成中...")
        
        # # 枚数決める
        # #println(data.t)
        # epoch_max = 100
        # epoch = len(self.sol.t)
        # if epoch < epoch_max
        #     step = 1
        # else
        #     step = div(epoch, epoch_max)
        
        # 軸揃える
        x_max = 0.2
        x_min = -0.2
        y_max = 0.2
        y_min = -0.2
        z_max = 0.55
        z_min = 0.0
        max_range = max(x_max-x_min, y_max-y_min, z_max-z_min)*0.5
        x_mid = (x_max + x_min) / 2
        y_mid = (y_max + y_min) / 2
        z_mid = (z_max + z_min) / 2

        self.xl = x_mid-max_range
        self.xu = x_mid+max_range
        self.yl = y_mid-max_range
        self.yu = y_mid+max_range
        self.zl = z_mid-max_range
        self.zu = z_mid+max_range
        
        
        self.time_template = 'time = %s [s]'
        
        
        fig = plt.figure()
        self.ax = fig.add_subplot(projection = '3d')
        
        ani = anm.FuncAnimation(
            fig = fig,
            func = self._update,
            frames = len(self.sol.t),
            #frames = int(self.TIME_SPAN / self.TIME_INTERVAL),
            interval = self.TIME_INTERVAL * 0.001
        )
        
        
        plt.show()

        ani.save("hoge.gif", fps=1/self.TIME_INTERVAL, writer='pillow')
        
        
        print("完了!")
        
        
        
        return None






if __name__ == "__main__":
    
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.plot()
    hoge.reproduce_state()
    hoge.make_aniation()