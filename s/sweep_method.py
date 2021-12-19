"""再現シミュレーション"""
from typing import Optional
import numpy as np
import matplotlib.pyplot as plt

import tqdm  # プログレスバーを出す


# パラメータの設定
LAMBDA = 271.23 # per day # 出生率
BETA_E0 = 3.11e-8 # person/day
BETA_I0 = 0.62e-8 # person/day
BETA_V0 = 1.03e-8 # person/day
C = 1.01e-4 # 係数
MU = 3.01e-5 # per day # 自然死亡率
ALPHA = 1/7 # 潜伏期間の逆数
RATE_W = 0.01 # per day #ウイルスの致死率
GAMMA = 1/15 # per day #回復率
SIGMA = 1 # per day # ウイルス除去率
XI_1 = 2.30 # 潜伏感染者からの排出率
XI_2 = 0 # per day # 感染者からの排出率
EPSILON = 0.05 # 制御入力の最大値調整係数

#　最適制御の重み係数
C1 = 0.25
C2 = 0.7
B1 = 8
B2 = 290 
B3 = 65



# 初期条件の定義
S_0 = 8998505
E_0 = 1000
I_0 = 475
R_0 = 10
V_0 = 10000
u1_0 = 0
u2_0 = 0
u3_0 = 0


# 区間の分割の設定
T = 300#19 # 開始日時(0日)から終端日時(19日) 
n = 1000#38 # 刻み幅 # 精度
step = T / n # 整数にしましょう
t = np.arange(0, T, step)  # 時間軸の値


def calc_beta(E, I, V):
    """非定常の遷移関数"""
    return (
        BETA_E0 / (1 + C * E),
        BETA_I0 / (1 + C * I),
        BETA_V0 / (1 + C * V),
    )


# 状態方程式 #引数と関数の数が一致するように
def f_ward(S,E,I,V,beta_E,beta_I,beta_V,u1,t=0):
    return LAMBDA - (1-u1) * (beta_E * S * E + beta_I * S * I + beta_V * S * V) - MU * S

def g_ward(S,E,I,V,beta_E,beta_I,beta_V,u1,t=0):
    return (1-u1) * (beta_E * S * E + beta_I * S * I + beta_V * S * V) - (ALPHA+MU)*E

def h_ward(E,I,u2,t=0):
    return ALPHA*E-(RATE_W+GAMMA+u2+MU)*I

def q_ward(I,R,u2,t=0):
    return (u2+GAMMA)*I-MU*R

def r_ward(E,I,V,u3,t=0):
    return XI_1*E-XI_2*I-(SIGMA+u3)*V



def one_step_of_state(X, U, step):
    """ルンゲクッタのワンステップ（状態）
    """
    
    def dX(X, U):
        S, E, I, R, V = X
        u1, u2, u3 = U
        beta_E,beta_I,beta_V = calc_beta(E, I, V)
        return np.array([
            f_ward(S,E,I,V,beta_E,beta_I,beta_V,u1),
            g_ward(S,E,I,V,beta_E,beta_I,beta_V,u1),
            h_ward(E,I,u2),
            q_ward(I,R,u2),
            r_ward(E,I,V,u3),
        ])

    k1 = dX(X, U)
    k2 = dX(X+k1*step/2, U)
    k3 = dX(X+k2*step/2, U)
    k4 = dX(X+k3*step, U)
    
    next_X = X + (k1 + 2*k2 + 2*k3 + k4) * step/6
    return next_X


# ハミルトニアンの微分
def f_back(E,I,V,u1,lam1,lam2,t=0):
    return (1 - u1) *\
        (BETA_E0 * E / (C * E + 1) + BETA_I0 * I / (C * I + 1) + BETA_V0 * V / ( C * V + 1)) *\
            (lam1-lam2) + lam1 * MU

def g_back(S,E,u1,lam1,lam2,lam3,lam5,t=0):
    return -C2 + (1-u1) *\
        (BETA_E0 * S / (C * E + 1) - BETA_E0 * C * S * E / (C*E+1)**2) *\
            (lam1-lam2) + lam2 * (ALPHA + MU) - lam3 * ALPHA - lam5 * XI_1

def h_back(S,I,u1,u2,lam1,lam2,lam3,lam4,lam5,t=0):
    return -C1 + (1 - u1) *\
        (BETA_I0 * S / (C * I + 1) - BETA_I0 * C * S * I / (C * I + 1)**2) * (lam1-lam2) +\
            lam3 * (RATE_W + GAMMA + MU + u2) -\
                lam4 * (u2 + GAMMA) - lam5 * XI_2

def q_back(lam4,t=0):
    return MU * lam4

def r_back(S,V,u1,u3,lam1,lam2,lam5,t=0):
    return (1-u1) * (BETA_V0 * S / ( C * V + 1) + BETA_V0 * C * S * V / (C * V + 1)**2) *\
        (lam1 - lam2) + lam5 * (SIGMA + u3)


def one_step_of_lambda(X, U, L, step):
    """逆方向ルンゲクッタのワンステップ（随伴ベクトル）"""
    
    def dL(L, X, U):
        lam1, lam2, lam3, lam4, lam5 = L
        S, E, I, R, V = X
        u1, u2, u3 = U
        
        z = np.array([
            f_back(E,I,V,u1,lam1,lam2),
            g_back(S,E,u1,lam1,lam2,lam3,lam5),
            h_back(S,I,u1,u2,lam1,lam2,lam3,lam4,lam5),
            q_back(lam4),
            r_back(S,V,u1,u3,lam1,lam2,lam5),
        ])
        return z
    
    k1 = dL(L, X, U)
    k2 = dL(L-k1*step/2, X, U)
    k3 = dL(L-k2*step/2, X, U)
    k4 = dL(L-k3*step, X, U)
    
    next_L = L - (k1 + 2*k2 + 2*k3 + k4) * step/6
    
    return next_L





# 結果を返すための配列(行列)の宣言 #初期条件を配列に追加
state = np.zeros((n, 5))
input = np.zeros((n, 3))
lagrange = np.zeros((n, 5))

# 初期値
state[0, :] = S_0, E_0, I_0, R_0, V_0
input[0, :] = u1_0, u2_0, u3_0
lagrange[-1, :] = 0, 0, 0, 0, 0


delta = 0.0001

saidai = 100  # whileだと無限ループする可能性がある．forのほうが良い気がする．

ts = []

for _i in tqdm.tqdm(range(saidai)):
    if _i > 0 and test > -1e-5:
        break

    old_state = state
    old_input = input
    old_lagrange = lagrange
    
    ### 順方向 ###
    for i in range(n-1):
        input_aver = (input[i, :] + input[i+1, :]) / 2
        state[i+1, :] = one_step_of_state(
            X = state[i, :],
            U = input_aver,
            step = step
        )

    ### 逆方向 ###
    for i in range(n-1):
        j = n -1 - i
        # 定義した関数の引数を全て持ってくること
        input_aver = (input[j, :] + input[j-1, :]) / 2
        lagrange[j-1, :] = one_step_of_lambda(
            X = state[j, :],
            U = input_aver,
            L = lagrange[j, :],
            step = step,
        )

    # 最適制御入力の配列 #_XによりX,lamの配列に要素が格納される
    # 配列の形で構成
    
    def calc_optiomal_u(X, L):
        """最適入力を計算"""
        optiomal_input = np.zeros((n, 3))
        for i in range(n):
            l1, l2, l3, l4, l5 = L[i, :]
            S, E, I, R, V = X[i, :]
            u1 = min(
                max(
                    0,
                    (l2-l1) / B1 *\
                        (BETA_E0*S*E/(C*E+1) + BETA_I0*S*I/(C*I+1) + BETA_V0*S*V/(C*V+1))
                ),
                1 - EPSILON,
            )
            u2 = min(max(0, (l3-l4)*I/B2), 1-EPSILON)
            u3 = min(max(0, l5*V/B3), 1-EPSILON)
            
            optiomal_input[i, :] = u1, u2, u3
        
        return optiomal_input
    
    optimal_input = calc_optiomal_u(state, lagrange)
    
    # 制御入力の調整 # lamの配列に要素が格納される # 配列同士の和の算出
    input = (optimal_input + old_input) / 2
    # c = 0.00000001
    # input = optimal_input * (1-c**_i) + old_input * c**_i

    # 収束条件の判定
    now = np.concatenate([state, input, lagrange], axis=1)
    old = np.concatenate([old_state, old_input, old_lagrange], axis=1)
    
    temp = delta * np.sum(np.abs(now), axis=0) - np.sum(np.abs(old - now), axis=0)
    
    #print(temp)
    test = min(temp)
    #print(test)
    ts.append(test)



# グラフで可視化
fig = plt.figure() # figureインスタンスを作成
ax1 = fig.add_subplot(221) # figureオブジェクトにaxesを追加#1行1列の1番目に導入する
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
ax1.plot(t, state[:, 1], color = "m", label = "Exposed population with control") # 潜伏感染者数をプロット
ax1.plot(t, state[:, 2], color = "k",label = "infected population with control") # 感染者数の曲線をプロット
ax2.plot(t, input[:, 0], color = "r", label = "u1") # 制御入力のプロット
ax2.plot(t, input[:, 1], color = "g",label = "u2")
ax2.plot(t, input[:, 2], color = "m",label = "u3")

ax4.plot(range(len(ts)), ts)
ax4.set_xlabel("traial")
ax4.set_ylabel("test score")
ax4.set_xlim(0, len(ts)-1)
ax4.set_ylim(min(ts), 0)
ax4.grid()

for i in range(5):
    ax3.plot(t, lagrange[:, i], label="lambda" + str(i+1))
ax3.set_xlabel("Time [days]")
ax3.set_ylabel('lagrange vector')
ax3.grid()
ax3.legend()
ax3.set_xlim(0,300)

# 実際の統計データをプロットする

ax1.set_xlabel('Time [days]') # x軸にラベルを追加
ax1.set_ylabel('Cumulative Comfirmes Case') # y軸にラベルを追加
ax2.set_xlabel('Time [days]') # x軸にラベルを追加
ax2.set_ylabel('control function') # y軸にラベルを追加
ax1.set_xlim(0,300) #　軸に範囲を追加
ax1.set_ylim(0,) #　軸に範囲を追加
ax2.set_xlim(0,300)
ax2.set_ylim(0, 1)

#ax.grid(True) # グリッドを入れる#網目の線
ax1.legend() # 凡例の追加
ax2.legend() # 凡例の追加
ax1.grid()
ax2.grid()
#ax.set_aspect('equal', adjustable='box')  # 軸を揃える#正方形 
plt.savefig('test.jpg',dpi=300)
plt.show()  # プロットを表示
