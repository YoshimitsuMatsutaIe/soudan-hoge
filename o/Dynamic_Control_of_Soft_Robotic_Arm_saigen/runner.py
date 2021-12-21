"""ゴミファイル"""


import numpy as np

from inverse_kinematics import InverseKinematics


hoge = InverseKinematics()

qd = hoge.calc_qd_from_xd(np.array([[0, 0.1, 0.147]]).T)
print(qd)