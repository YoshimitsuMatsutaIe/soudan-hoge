"""使用例


"""

import numpy as np


from simulator import Simulator
from inverse_kinematics import InverseKinematics


def pd(t):
    """タスク空間上の所望の位置"""
    return np.array([
        [0.1 * np.sin(3*t)],
        [0.1 * np.cos(3*t)],
        [0.147],
    ])


def pd_dot(t):
    """タスク空間上の所望の測度"""
    return np.array([
        [0.1 * 3 * np.cos(3*t)],
        [0.1 * 3 * -np.sin(3*t)],
        [0],
    ])


def pd_dot_dot(t):
    """タスク空間上の所望の加速度"""
    return np.array([
        [0.1 * 9 * -np.sin(3*t)],
        [0.1 * 9 * -np.cos(3*t)],
        [0],
    ])



# hoge = Simulator(
#     TIME_SPAN = 2*np.pi/3,
#     TIME_INTERVAL = 0.01,
#     pd = pd,
#     pd_dot = pd_dot,
#     pd_dot_dot = pd_dot_dot
# )

# hoge.run_simulation()
# #hoge.plot_actuator_data()
# hoge.plot_all()
# hoge.make_animation()


hoge2 = Simulator()
hoge2.plot_test(np.array([[0.011, 0.01, 0.01]]).T)



# wow = InverseKinematics()

# a, b, c = wow.calc_desired_actuator_state_all(
#     pd, pd_dot, pd_dot_dot, 1, 0.01
# )


# import matplotlib.pyplot as plt
# t = np.arange(0, 1, 0.01)
# fig = plt.figure()
# ax = fig.add_subplot(1, 3, 1)
# ax.plot(t, a[:, 0], label="l1_d")
# ax.plot(t, a[:, 1], label="l2_d")
# ax.plot(t, a[:, 2], label="l3_d")
# ax.legend()
# ax.grid()

# ax1 = fig.add_subplot(1, 3, 2)
# ax1.plot(t, b[:, 0], label="l1_dot_d")
# ax1.plot(t, b[:, 1], label="l2_dot_d")
# ax1.plot(t, b[:, 2], label="l3_dot_d")
# ax1.legend()
# ax1.grid()

# ax2 = fig.add_subplot(1, 3, 3)
# ax2.plot(t, c[:, 0], label="l1_dot_dot_d")
# ax2.plot(t, c[:, 1], label="l2_dot_dot_d")
# ax2.plot(t, c[:, 2], label="l3_dot_dot_d")
# ax2.legend()
# ax2.grid()

# plt.show()