import ctypes as ct
import numpy as np


#/home/matsuta/src/ctrlab2021_soudan/o/soft_robot/derivation_of_kinematics/derived/so
f = ct.cdll.LoadLibrary("o/soft_robot/derivation_of_kinematics/derived/so/Phi_0.so")
ff = f.Phi_0
ff.argtypes = (ct.POINTER(ct.c_longdouble), ct.c_longdouble, ct.POINTER(ct.c_longdouble))
ff.restype= ct.POINTER(ct.c_double)

q = np.zeros((3, 1)).astype(np.float64)
x = np.zeros((3, 1)).astype(np.float64)

qc = q.ctypes.data_as(ct.POINTER((ct.c_longdouble * 3) * 1)).contents
xc = x.ctypes.data_as(ct.POINTER((ct.c_longdouble * 3) * 1)).contents

ff(q, 1.0, x)
print(Y)