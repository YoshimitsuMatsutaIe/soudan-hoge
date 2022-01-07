import numpy
def f(q, q_dot, xi):
    l1, l2, l3 = q[0,0], q[1,0], q[2,0]
    l1_dot, l2_dot, l3_dot = q_dot[0,0], q_dot[1,0], q_dot[2,0]

    return numpy.array([[1700, 0, 0], [0, 1700, 0], [0, 0, 1700]])