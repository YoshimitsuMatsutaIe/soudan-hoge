"""シミュレーション関係

上から順にメソッドを実行してけば良い  

"""

import numpy as np
from math import pi, cos, sin
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as anm
from matplotlib.patches import Circle
import mpl_toolkits.mplot3d.art3d as art3d
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
import yaml

from kinematics import Kinematics
from differential_kinematics import DifferentialKinematics
from controller import OriginalRMPAttractor, OriginalRMPJointLimitAvoidance, pullback, PDFeedBack


def rotate_3d(alpha, beta, gamma):
    """3次元回転行列"""
    
    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(alpha), -np.sin(alpha)],
        [0, np.sin(alpha), np.cos(alpha)],
    ])
    Ry = np.array([
        [np.cos(beta), 0, np.sin(beta)],
        [0, 1, 0],
        [-np.sin(beta), 0, np.cos(beta)],
    ])
    Rz = np.array([
        [np.cos(gamma), -np.sin(gamma), 0],
        [np.sin(gamma), np.cos(gamma), 0],
        [0, 0, 1],
    ])
    
    return Rx @ Ry @ Rz



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
    
    ### デフォルト値 ###
    TIME_SPAN = 10  # シミュレーション時間
    TIME_INTERVAL = 0.01  # 刻み時間
    
    env_param = {
        "goal_param" : {4 : lambda t: np.array([[0.2, 0.1, 0.6]])},
        "goal_dot_param" : None,
        "goal_dot_dot_param" : None,
        "target_param" : None,
    }
    
    controller_param = {
        "name" : "rmp",
        "attractor" : {
            "max_speed" : 10,
            "gain" : 300,
            "a_damp_r" : 0.05,
            "sigma_W" : 1,
            "sigma_H" : 1,
            "A_damp_r" : 0.01,
        },
        "jlavoidance" : {
            "gamma_p" : 5,
            "gamma_d" : 1,
            "lam" : 1
        }
    }
    
    
    def __init__(self, N=5, env_param=None, controller_param=None,):
        """
        Parameters
        ---
        N : int
            セクションの数
        env_param : dict[int, list[float]]
            環境のパラメータ．キーが対象となるセクション，値が目標位置ベクトル
        controller_param : dict[str, float]
            コントローラのパラメータ
        """
        
        self.N = N  # セクションの数
        
        self.kim = Kinematics()
        self.diff_kim = DifferentialKinematics(N)
        
        
        if env_param is not None:
            self.env_param = env_param
        self.set_environment()
        
        if controller_param is not None:
            self.controller_param = controller_param
        self.set_controller()


    def set_environment(self,):
        """シミュレーション環境をセット"""
        
        self.goal = self.env_param["goal_param"]
        self.goal_dot = self.env_param["goal_dot_param"]
        self.goal_dot_dot = self.env_param["goal_dot_dot_param"]


    def set_controller(self,):
        """制御器をセット"""
        
        if self.controller_param["name"] == "pdfb":
            self.pdfb = PDFeedBack(**self.controller_param)
            self.calc_input = self.calc_input_by_PDFB
        
        else self.controller_param["name"] == "rmp":
            self.q_max = np.array([[0.1] * self.N*3]).T
            self.q_min = -self.q_max
            
            self.attractor = OriginalRMPAttractor(
                **self.controller_param["attractor"]
            )
            self.jlavoidance = OriginalRMPJointLimitAvoidance(
                **self.controller_param["jlavoidance"]
            )
            self.calc_input = self.calc_input_by_RMP

    
    
    def calc_input_by_PDFB(self, t, q, q_dot):
        """RMPで制御入力（角加速度）を計算"""
        
        return self.pdfb.input(
            x = self.kim.Phi(self.N, q),
            J = self.diff_kim.J(self.N, q),
            x_dot = J * q_dot,
            xd = self.goal[self.N](t)
        )
    
    
    def calc_input_by_RMP(self, t, q, q_dot):
        """RMPで制御入力（角加速度）を計算"""
        
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
                xs[k], x_dots[k], self.goal[k](t)
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
        return np.linalg.pinv(root_M) @ root_f
    


    def X_dot(self, t, X):
        """scipyで使う微分方程式"""
        
        #print("t = ", t)
        q = X[:self.N*3].reshape(self.N*3, 1)
        q_dot = X[self.N*3:].reshape(self.N*3, 1)
        
        return np.ravel(
            np.concatenate(
                [
                    q_dot,
                    self.calc_input(t, q, q_dot)
                ]
            )
        )
    
    
    @stop_watch
    def run(self, TIME_SPAN=None, TIME_INTERVAL=None):
        """シミュレーション実行"""
        
        date_now = datetime.datetime.now()  # 名前つけるとき使う
        name = date_now.strftime('%Y-%m-%d--%H-%M-%S')
        cwd = str(Path().resolve())
        self.base = cwd + "/result/" + name
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

        print("xd...")
        self.xd_data = []
        for i in tqdm.tqdm(range(len(self.sol.t))):
            temp = []
            for j in range(self.N):
                temp.append(self.goal)
            self.q_dot_data.append(np.array([temp]).T)


    @stop_watch
    def save_data(self,):
        """いろいろデータ保存"""
        
        print("データ保存中...")
        
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

        # シミュレーション設定
        with open(self.base + "/envireonment.yaml", "w") as f:
            yaml.dump(self.env_param, f)
        
        # 使ったパラメータ保存
        with open(self.base + "/controller_param.yaml", "w") as f:
            yaml.dump(self.controller_param, f)
        
        # with open(self.base + "/jlavoidance_param.yaml", "w") as f:
        #     yaml.dump(self.jlavoidance_param, f)
        
        
        print("完了!")



    @stop_watch
    def plot_basic(self,):
        """基本的なものをプロット"""
        
        fig = plt.figure(figsize=(14, 15))
        
        # 誤差の時系列
        es = []
        for k in self.goal:
            e = []
            for t, i in enumerate(self.sol.t):
                e.append(
                    np.linalg.norm(self.goal[k](t) - self.x_data[i][k][:, -1])
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






    
    def _update(self, i):
        """アニメのフレーム作成"""
        
        t = i * self.TIME_INTERVAL
        
        
        self.ax.cla()
        
        self.ax.grid(True)
        self.ax.set_xlabel('X [m]')
        self.ax.set_ylabel('Y [m]')
        self.ax.set_zlabel('Z [m]')
        
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
                s = 200, label = 'Goal of Sec ' + str(k),
                marker = '*',# color = '#ff7f00', 
                alpha = 1, linewidths = 1.1, edgecolors = 'red'
            )
        
        # アーム
        for j in range(self.N):
            self.ax.plot(
                self.x_data[i][j][0,:], self.x_data[i][j][1,:], self.x_data[i][j][2,:],
                label="Section" + str(j)
            )

        # 接続位置
        for j in range(self.N):
            self.ax.scatter(
                self.x_data[i][j][0,-1], self.x_data[i][j][1,-1], self.x_data[i][j][2,-1],
            )
        # p = Circle(
        #     xy = (0.1, 0.1),
        #     radius = 0.0125,
        #     angle=0.4,
        #     alpha = 0.5 
        # )
        # self.ax.add_patch(p)
        # art3d.pathpatch_2d_to_3d(p, z=0, zdir = 'x')
        #art3d.pathpatch_translate(p, (0.5, 1, 0))

        Xc,Yc,Zc = self.data_for_cylinder_along_z(
            0.2, 0, 0.05, 0.7
        )
        self.ax.plot_surface(Xc, Yc, Zc, alpha=0.4)

        self.ax.legend()


    def data_for_cylinder_along_z(self, center_x,center_y,radius,height_z):
        """コピペ"""
        z = np.linspace(0, height_z, 50)
        theta = np.linspace(0, 2*np.pi, 50)
        theta_grid, z_grid=np.meshgrid(theta, z)
        x_grid = radius*np.cos(theta_grid) + center_x
        y_grid = radius*np.sin(theta_grid) + center_y
        return x_grid,y_grid,z_grid


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
        # ani.save(
        #     self.base + "/animation.mp4", fps=60, writer='ffmpeg'
        # )  #Windowsだと無理かも
        
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


# 就活
def ex_yasukawa():
    """安川電機プレゼン用
    
    ・研究概要 動画不可
    """
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.reproduce_state()
    hoge.plot_basic()
    hoge.make_aniation()


def ex_denso():
    """デンソー用
    
    ・スライド．動画可?
    """
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.reproduce_state()
    hoge.plot_basic()
    hoge.make_aniation()


def ex_kawasaki():
    """川崎重工プレゼン用
    
    
    """
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.reproduce_state()
    hoge.plot_basic()
    hoge.make_aniation()


def ex_fujikoshi():
    """不二越 研究プレゼン
    
    動画不可？
    """
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.reproduce_state()
    hoge.plot_basic()
    hoge.make_aniation()



if __name__ == "__main__":
    
    hoge = Simulator(N=5, goal_param=None)
    
    hoge.run(TIME_SPAN = 10)
    hoge.reproduce_state()
    hoge.save_data()
    hoge.plot_basic()
    hoge.make_aniation()