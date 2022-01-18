import numpy as np

from math import sin, cos, sqrt

import matplotlib.pyplot as plt

import sympy
import scipy
import scipy.optimize
from scipy.optimize import fsolve
from scipy.integrate import odeint


TIME_SPAN = 5.2
dt = 0.1


T = np.arange(0.0, 5.2, 0.1)
T1 = np.arange(0.0, 5.1, 0.1)
T2 = np.arange(0.0, 5.0, 0.1)


p1 = 2106.99588477366
p2 = 1053.49794238683
p3 = 948.148148148148
p4 = 6320.98765432099
p5 = 1422.22222222222
p6 = 8.88888888888889
p7 = 158.024691358025
p8 = 3160.49382716049
p9 = 4.44444444444444
p10 = 474.074074074074
p11 = 71.1111111111111
p12 = 0.00299999999999997
p13 = sqrt(3)



def qd(L, t):
    """[l1, l2, l3]から参照軌道qdを求める"""
    
    l1, l2, l3 = L
    
    qd_1 = p1*l1**4 - p2*l1**3*l2 - p2*l1**3*l3 + p3*l1**3 - p4*l1**2*l2*l3  - \
        p5*l1**2*l2 - p5*l1**2*l3 - p6*l1**2 + p1*l1*l2**3 + p8*l1*l2**2*l3  + \
            p5*l1*l2**2 + p8*l1*l2*l3**2 - p9*l1*l2 + p1*l1*l3**3 + p5*l1*l3**2 - \
                p9*l1*l3 - 4.0*l1 - p2*l2**4 - p2*l2**3*l3 - p10*l2**3 + p9*l2**2 -\
                        p2*l2*l3**3 + p6*l2*l3 + 2.0*l2 - p2*l3**4 - p10*l3**3 + p9*l3**2 + \
                            2.0*l3 - 0.1*sin(3*t)
    qd_2 = p2*p13*l1**3*l2 - p2*p13*l1**3*l3 + p10*p13*l1**2*l2 - p10*p13*l1**2*l3  -\
            p8*p13*l1*l2**2*l3 - p10*p13*l1*l2**2 + p8*p13*l1*l2*l3**2 - p9*p13*l1*l2 + \
                p10*p13*l1*l3**2 + p9*p13*l1*l3 + p2*p13*l2**4 - p2*p13*l2**3*l3 + \
                    p10*p13*l2**3  - p3*p13*l2**2*l3 - p9*p13*l2**2 + p2*p13*l2*l3**3 + \
                        p3*p13*l2*l3**2 - 2.0*p13*l2 - p2*p13*l3**4 - p10*p13*l3**3 + \
                            p9*p13*l3**2 + 2.0*p13*l3 - 0.1*cos(3*t)
    qd_3 = -p7*l1**3 - p11*l1**2 + p10*l1*l2*l3 + p11*l1*l2 + p11*l1*l3 + l1/3 - \
        p7*l2**3 - p11*l2**2 + p11*l2*l3 + l2/3 - p7*l3**3 - p11*l3**2 + l3/3 + p12
    
    return np.array([qd_1, qd_2, qd_3])


# #qdの微分
# def qd_dot(p, t):
#     return (fsolve(qd,[0, 0, 0],args=(t+0.1,))-fsolve(qd,[0, 0, 0],args=(t,)))/0.1



# #qdの2階微分
# def qd_dot_dot(p, t):
#     return (((fsolve(qd,[0, 0, 0],args=(t+0.2,))-fsolve(qd,[0, 0, 0],args=(t+0.1,)))/0.1)-((fsolve(qd,[0, 0, 0],args=(t+0.1,))-fsolve(qd,[0, 0, 0],args=(t,)))/0.1))/0.1



#qの算出
def X_dot(X, t):
    """微分方程式"""
    
    L = X[:3]
    L_dot = X[3:]
    
    
    dl0=(fsolve(qd,[0, 0, 0],args=(t+0.1,))-fsolve(qd,[0, 0, 0],args=(t,)))/0.1
    
    d2l0= (((fsolve(qd,[0, 0, 0],args=(t+0.2,))-fsolve(qd,[0, 0, 0],args=(t+0.1,)))/0.1)-((fsolve(qd,[0, 0, 0],args=(t+0.1,))-fsolve(qd,[0, 0, 0],args=(t,)))/0.1))/0.1
    
    #print("dl0 = ", dl0)
    #print("d2l0 = ", d2l0)
    
    dq0dt = np.concatenate(
        [
            L_dot,
            -200*L_dot - 10000*L + 200*dl0 + 10000*d2l0
        ],
        axis = 0
    )
    #print("dq0dt = ", dq0dt)
    return dq0dt



X0 = np.zeros(6)  # 初期値


so12 = odeint(
    func = X_dot,
    y0 = X0,
    t=T2
) #微分方程式を解く


# so12 = odeint(
#     fun,
#     np.zeros(6),
#     T2
# ) #微分方程式を解く

fig = plt.figure()
ax = fig.add_subplot()
ax.plot(so12)



plt.show()