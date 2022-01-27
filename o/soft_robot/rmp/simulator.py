"""シミュレーション関係

上から順にメソッドを実行してけば良い  

"""

import numpy as np
from math import pi, cos, sin
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import matplotlib.animation as anm
import matplotlib.patches as patches
import mpl_toolkits.mplot3d.art3d as art3d
#from matplotlib.font_manager import FontProperties
#plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.family"] = "DejaVu Serif"
from IPython.display import HTML
from tqdm import tqdm
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
from goal_trajectory import Point, Circle, RoseCurve


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
    TIME_SPAN = 3  # シミュレーション時間
    TIME_INTERVAL = 0.01  # 刻み時間
    
    #dddd =  np.array([[0.2, 0.1, 0.6]]).
    env_param = {
        "goal" : {
            4 : {
                "name" : "Circle",
                "param" :{
                    "r" : 0.2,
                    'center' : [0.0, 0.0, 0.75],
                    "omega" : 1.5,
                    "alpha" : 0,
                    "beta" : 0,
                    "gamma" : 0
                }
            }
        },
        "target_param" : None,
    }
    
    controller_param = {
        "name" : "rmp",
        "attractor" : {
            "max_speed" : 1800,
            "gain" : 15000,
            "a_damp_r" : 0.05,
            "sigma_W" : 1,
            "sigma_H" : 1,
            "A_damp_r" : 0.01,
        },
        "jlavoidance" : {
            "gamma_p" : 4,
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
            環境のパラメータ
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
        
        # 目標位置
        self.goal = {}
        
        for k, v in self.env_param["goal"].items():
            if v['name'] == 'Point':
                self.goal[k] = Point(**v['param'])
            elif v['name'] == 'Circle':
                self.goal[k] = Circle(**v['param'])
            elif v['name'] == 'RoseCurve':
                self.goal[k] == RoseCurve(**v['param'])


    def set_controller(self,):
        """制御器をセット"""
        
        if self.controller_param["name"] == "pdfb":
            print("PDフィードバック制御")
            self.pdfb = PDFeedBack(**self.controller_param)
            self.calc_input = self.calc_input_by_PDFB
        
        elif self.controller_param["name"] == "rmp":
            print("RMP制御")
            
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
            xd = self.goal[self.N].xd(t)
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
        for k, v in self.goal.items():
            f, M = self.attractor.get_natural(
                x = xs[k],
                dx = x_dots[k],
                x0 = v.xd(t),
                dx0 = v.xd_dot(t),
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
        
        print("t = ", t)
        q = X[:self.N*3].reshape(self.N*3, 1)
        q_dot = X[self.N*3:].reshape(self.N*3, 1)
        
        return np.ravel(np.concatenate(
            [
                q_dot,
                self.calc_input(t, q, q_dot)
            ]
        ))
    
    
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
        for i in tqdm(range(len(self.sol.t))):
            temp = []
            for j in range(self.N):
                temp.append(self.sol.y[3*j][i])
                temp.append(self.sol.y[3*j+1][i])
                temp.append(self.sol.y[3*j+2][i])
            self.q_data.append(np.array([temp]).T)
        
        print("q_dot...")
        self.q_dot_data = []
        for i in tqdm(range(len(self.sol.t))):
            temp = []
            for j in range(self.N):
                temp.append(self.sol.y[3*j+3*self.N][i])
                temp.append(self.sol.y[3*j+1+3*self.N][i])
                temp.append(self.sol.y[3*j+2+3*self.N][i])
            self.q_dot_data.append(np.array([temp]).T)
        
        print("x...")
        self.x_data = []
        for i in tqdm(range(len(self.sol.t))):
            temp = []
            for j in range(self.N):
                temp.append(
                    np.concatenate(
                        [self.kim.Phi(j+1, self.q_data[i], xi) for xi in self.kim.xi_large],
                        axis=1
                    )
                )
            self.x_data.append(temp)

        print("xd and error...")
        self.xd_data = []
        self.error_data = []
        for i, t in enumerate(tqdm(self.sol.t)):
            temp_xd = []
            temp_error = []
            for k, v in self.goal.items():
                xd = v.xd(t)
                temp_xd.append(xd)
                temp_error.append(np.linalg.norm(xd - self.x_data[i][k][:, [-1]]))
            self.xd_data.append(temp_xd)
            self.error_data.append(temp_error)


    @stop_watch
    def save_data(self,):
        """いろいろデータ保存
        
        もっと賢いやり方があると思います
        """
        
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
        
        self.actuator_data = np.concatenate(_temp, axis=1)
        np.savetxt(
            self.base + '/actuator.csv',
            self.actuator_data,
            header=header,
            comments='',
            delimiter = ","
        )
        
        # x（エンドエフェクタのみ）, xd, errorを保存
        n = ('x', 'y', 'z')
        header = 't'
        for i in range(self.N):
            for j in range(3):
                header += ',p_' + str(i) + '_' + n[j]
        for i in self.goal.keys():
            for j in range(3):
                header += ',pd_' + str(i) + '_' + n[j]
        for i in self.goal.keys():
            header += ',error_of_sec_' + str(i)
        
        temp = []
        for i, t in enumerate(self.sol.t):
            _temp = [t.reshape(1, 1)]
            for j in range(self.N):
                _temp.append(self.x_data[i][j][:, [-1]])
            for j in range(len(self.goal)):
                _temp.append(self.xd_data[i][j])
            for j in range(len(self.goal)):
                _temp.append(np.array([[self.error_data[i][j]]]))
            temp.append(np.concatenate(_temp))
        
        temp = np.concatenate(temp, axis=1)
        self.task_data = temp.T
        np.savetxt(
            self.base + '/task.csv',
            self.task_data,
            header=header,
            comments='',
            delimiter = ","
        )

        # シミュレーション設定を保存
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
        
        ax = fig.add_subplot(3, 1, 1)
        _max = []
        for i, k in enumerate(self.goal.keys()):
            ax.plot(
                self.sol.t, self.task_data[:, 1 + 3*self.N + 3*len(self.goal) + i],
                label=r"Section " + str(k)
            )
            _max.append(np.max(self.task_data[:, 1 + 3*self.N + 3*len(self.goal) + i]))
        ax.set_xlabel(r'Time $\it{t}$ [s]')
        ax.set_ylabel(r"Position Error to Goal [m]")
        ax.set_xlim(0, self.TIME_SPAN)
        ax.set_ylim(0, max(_max))
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
        return



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
        for k, v in self.goal.items():
            xd = v.xd(t)
            self.ax.scatter(
                xd[0, 0], xd[1, 0], xd[2, 0],
                s = 200, label = 'Goal of Sec ' + str(k),
                marker = '*',# color = '#ff7f00', 
                alpha = 1, linewidths = 1.1, edgecolors = 'red'
            )
        
        
        # 目標の軌道
        for j, k in enumerate(self.goal.keys()):
            _i0 = 1 + 3*self.N + 3*j
            self.ax.plot(
                self.task_data[:, _i0],
                self.task_data[:, _i0+1],
                self.task_data[:, _i0+2],
                label = 'Goal traj of Sec ' + str(k),
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

        # Xc,Yc,Zc = self.data_for_cylinder_along_z(
        #     0.2, 0, 0.05, 0.7
        # )
        # self.ax.plot_surface(Xc, Yc, Zc, alpha=0.4)

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
        
        # 枚数決める
        #println(data.t)
        epoch_max = 100
        epoch = len(self.sol.t)
        if epoch < epoch_max:
            step = 1
        else:
            step = epoch // epoch_max
        
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
            #frames = epoch_max,
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



def en_tuiju():
    """円軌道追従
    """
    
    env_param = {
        "goal_param" : {4 : lambda t: np.array([[0.2, 0.1, 0.6]]).T},
        "goal_dot_param" : None,
        "goal_dot_dot_param" : None,
        "target_param" : None,
    }
    
    sim = Simulator(N=5)
    



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
    
    hoge = Simulator(N=5)
    
    hoge.run(TIME_SPAN = 5)
    hoge.reproduce_state()
    hoge.save_data()
    hoge.plot_basic()
    hoge.make_aniation()