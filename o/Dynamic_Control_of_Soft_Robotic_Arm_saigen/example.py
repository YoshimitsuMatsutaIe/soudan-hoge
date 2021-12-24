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
#     TIME_SPAN = 2,
#     TIME_INTERVAL = 0.01,
#     pd = pd,
#     pd_dot = pd_dot,
#     pd_dot_dot = pd_dot_dot
# )

# hoge.run_simulation()
# #hoge.plot_actuator_data()
# hoge.plot_all()
# #hoge.make_animation()


# hoge2 = Simulator()
# hoge2.plot_test(np.array([[0.01, 0, 0]]).T)



wow = InverseKinematics()

a, b, c = wow.calc_desired_actuator_state_all(
    pd, pd_dot, pd_dot_dot, 1, 0.01
)
print(a)