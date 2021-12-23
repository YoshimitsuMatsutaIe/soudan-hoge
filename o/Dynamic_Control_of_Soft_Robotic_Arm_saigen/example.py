"""使用例


"""

import numpy as np
from math import sin, cos, sqrt, pi


from simulator import Simulator



def pd(t):
    """タスク空間上の所望の位置"""
    return np.array([
        [0.1 * sin(3*t)],
        [0.1 * cos(3*t)],
        [0.147],
    ])


def pd_dot(t):
    """タスク空間上の所望の測度"""
    return np.array([
        [0.1 * 3 * cos(3*t)],
        [0.1 * 3 * -sin(3*t)],
        [0],
    ])


def pd_dot_dot(t):
    """タスク空間上の所望の加速度"""
    return np.array([
        [0.1 * 9 * -sin(3*t)],
        [0.1 * 9 * -cos(3*t)],
        [0],
    ])



hoge = Simulator(
    TIME_SPAN = 1,
    TIME_INTERVAL = 0.01,
    pd = pd,
    pd_dot = pd_dot,
    pd_dot_dot = pd_dot_dot
)

hoge.run_simulation()
hoge.plot_actuator_data()
hoge.make_animation()


# hoge2 = Simulator()
# hoge2.plot_test(np.array([[0.01, 0, 0]]).T)