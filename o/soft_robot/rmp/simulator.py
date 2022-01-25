import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as anm
#from matplotlib.font_manager import FontProperties
#plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.family"] = "DejaVu Serif"
from IPython.display import HTML
import tqdm
import time
import datetime
from functools import wraps
import os
from pathlib import Path
import csv

from kinematics import Kinematics
from differential_kinematics import DifferentialKinematics
from controller import OriginalRMPAttractor, OriginalRMPJointLimitAvoidance, pullback


def stop_watch(func):
    """時間計測"""
    @wraps(func)
    def wrapper(*args, **kargs):
        start = time.time()
        result = func(*args, **kargs)
        elapsed_time = time.time() - start

        print("{} s in {}".format(elapsed_time, func.__name__))
        return result
    return wrapper


class Simulator:
    """シミュレーション関係"""
    
    TIME_SPAN = 10  # デフォルト
    TIME_INTERVAL = 0.01  # デフォルト
    
    
    def __init__(self, N=3, goal=None):
        
        self.N = N  # セクションの数
        
        self.kim = Kinematics()
        self.diff_kim = DifferentialKinematics(N)
        
        
        # 目標位置
        if goal is None:  # デフォルト値
            self.goal = {
                4:np.array([[0.1, 0.1, 0.6]]).T,
                3:np.array([[0.1, 0.1, 0.5]]).T,
                2:np.array([[0.1, 0.1, 0.4]]).T,
            }
        else:
            self.set_goal(goal)
        
        # アクチュエータ制約
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
        
        

    def set_goal(self, goal):
        """目標位置を準備
        
        goal : dict
        """
        self.goal = goal
    
    
    def X_dot(self, t, X):
        """scipyで使う微分方程式"""
        #print("t = ", t)
        q = X[:self.N*3].reshape(self.N*3, 1)
        q_dot = X[self.N*3:].reshape(self.N*3, 1)
        
        # pushforward演算
        xs = [self.kim.Phi(i+1, q) for i in range(self.N)]
        Js = [self.diff_kim.J(i+1, q) for i in range(self.N)]
        x_dots = [J @ q_dot for J in Js]

        
        # ルートのRMP
        root_f = np.zeros((self.N*3, 1))
        root_M = np.zeros((self.N*3, self.N*3))


        # アトラクター
        for k in self.goal:
            f, M = self.attractor.get_natural(
                xs[k], x_dots[k], self.goal[k], np.zeros((3,1))
            )
            pullbacked_f, pullbacked_M = pullback(f, M, Js[k])
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
    
    
    @stop_watch
    def run(self, TIME_SPAN=None, TIME_INTERVAL=None):
        """シミュレーション実行"""
        
        date_now = datetime.datetime.now()  # 名前つけるとき使う
        name = date_now.strftime('%Y-%m-%d--%H-%M-%S')
        cwd = str(Path().resolve())
        self.base = cwd + "/" + name
        os.makedirs(self.base, exist_ok=True)
        
        
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
    
    
    @stop_watch
    def plot_basic(self,):
        """基本的なものをプロット"""
        
        fig = plt.figure(figsize=(14, 15))
        
        # 誤差の時系列
        es = []
        for k in self.goal:
            e = []
            for i in range(len(self.sol.t)):
                e.append(
                    np.linalg.norm(self.goal[k] - self.x_data[i][k][:, -1])
                )
            es.append(e)
        
        ax = fig.add_subplot(3, 1, 1)
        for i, k in enumerate(self.goal):
            ax.plot(self.sol.t, es[i], label=r"Section " + str(k))
        ax.set_xlabel(r'Time $\it{t}$ [s]')
        ax.set_ylabel(r"Position Error to Goal [m]")
        ax.set_xlim(0, self.TIME_SPAN)
        ax.set_ylim(0, max([x for row in es for x in row]))
        ax.legend()
        ax.grid()
        
        ax2 = fig.add_subplot(3, 1, 2)
        for i in range(self.N):
            ax2.plot(self.sol.t, self.sol.y[3*i], label=r"$\it{l^" + str(i) + r"_1}$")
            ax2.plot(self.sol.t, self.sol.y[3*i+1], label=r"$\it{l^" + str(i) + r"_2}$")
            ax2.plot(self.sol.t, self.sol.y[3*i+2], label=r"$\it{l^" + str(i) + r"_3}$")
        ax2.set_xlabel(r'Time $\it{t}$ [s]')
        ax2.set_ylabel(r'Actuator Length $\it{l}$ [m]')
        ax2.set_xlim(0, self.TIME_SPAN)
        ax2.legend()
        ax2.grid()
        
        ax3 = fig.add_subplot(3, 1, 3)
        for i in range(self.N):
            ax3.plot(self.sol.t, self.sol.y[3*i+3*self.N], label=r"$\.{l^" + str(i) + r"_1}$")
            ax3.plot(self.sol.t, self.sol.y[3*i+1+3*self.N], label=r"$\.{l^" + str(i) + r"_2}$")
            ax3.plot(self.sol.t, self.sol.y[3*i+2+3*self.N], label=r"$\.{l^" + str(i) + r"_3}$")
        
        ax3.set_xlabel(r'Time $\it{t}$ [s]')
        ax3.set_ylabel(r"Actuator Velocity $\.{l}$ [m/s]")
        ax3.set_xlim(0, self.TIME_SPAN)
        ax3.legend()
        ax3.grid()
        
        fig.savefig(self.base + "/basic.png")
        
        
        #plt.show()
        
        print("plot完了!")


    @stop_watch
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
                        [self.kim.Phi(j+1, self.q_data[i], xi) for xi in self.kim.xi_large],
                        axis=1
                    )
                )
            self.x_data.append(temp)


    @stop_watch
    def save_data(self,):
        """csvにデータを保存"""
        
        print("csvに保存中...")
        
        # q, q_dotを保存
        header = 't'
        for i in range(self.N):
            for j in range(3):
                header += ',l_' + str(i) + '_' + str(j)
        for i in range(self.N):
            for j in range(3):
                header += ',l_dot_' + str(i) + '_' + str(j)
        
        _temp = [self.sol.t.reshape(len(self.sol.t), 1)]
        for i in range(self.N*3*2):
            _temp.append(self.sol.y[i].reshape(len(self.sol.t), 1))
        
        _temp = np.concatenate(_temp, axis=1)
        np.savetxt(
            self.base + '/actuator.csv',
            _temp,
            header=header,
            comments='',
            delimiter = ","
        )
        
        print("完了!")


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
        for k in self.goal:
            self.ax.scatter(
                self.goal[k][0, 0], self.goal[k][1, 0], self.goal[k][2, 0],
                s = 200, label = 'goal of sec' + str(k),
                marker = '*',# color = '#ff7f00', 
                alpha = 1, linewidths = 1.1, edgecolors = 'red'
            )
        
        # アーム
        for j in range(self.N):
            self.ax.plot(
                self.x_data[i][j][0,:], self.x_data[i][j][1,:], self.x_data[i][j][2,:],
                label="section" + str(j)
            )

        # 接続位置
        for j in range(self.N):
            self.ax.scatter(
                self.x_data[i][j][0,-1], self.x_data[i][j][1,-1], self.x_data[i][j][2,-1],
            )

        self.ax.legend()



    @stop_watch
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
        z_max = 0.8
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
        
        
        fig = plt.figure(figsize=(8, 8))
        self.ax = fig.add_subplot(projection = '3d')
        
        ani = anm.FuncAnimation(
            fig = fig,
            func = self._update,
            frames = len(self.sol.t),
            #frames = int(self.TIME_SPAN / self.TIME_INTERVAL),
            interval = self.TIME_INTERVAL * 0.001
        )
        
        
        # 保存
        ani.save(
            self.base + "/animation.gif", fps=60, writer='pillow'
        )
        ani.save(
            self.base + "/animation.mp4", fps=60, writer='ffmpeg'
        )  #Windowsだと無理かも
        
        print("完了!")
        plt.show()
        
        return ani



def ex_default():
    """デフォルトの実行例"""
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.reproduce_state()
    hoge.plot_basic()
    hoge.make_aniation()





if __name__ == "__main__":
    
    hoge = Simulator(N=5)
    
    hoge.run(TIME_SPAN = 0.03)
    hoge.reproduce_state()
    hoge.save_data()
    hoge.plot_basic()
    hoge.make_aniation()