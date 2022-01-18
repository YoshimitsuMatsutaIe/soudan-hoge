import sympy as sy


t = sy.Symbol("t")
l1, l2, l3 = sy.symbols("l1, l2, l3")
p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13 = sy.symbols(
    "p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13"
)
sin_3t, cos_3t = sy.symbols("sin_3t, cos_3t")


qd_1 = p1*l1**4 - p2*l1**3*l2 - p2*l1**3*l3 + p3*l1**3 - p4*l1**2*l2*l3  - \
    p5*l1**2*l2 - p5*l1**2*l3 - p6*l1**2 + p1*l1*l2**3 + p8*l1*l2**2*l3  + \
        p5*l1*l2**2 + p8*l1*l2*l3**2 - p9*l1*l2 + p1*l1*l3**3 + p5*l1*l3**2 - \
            p9*l1*l3 - 4.0*l1 - p2*l2**4 - p2*l2**3*l3 - p10*l2**3 + p9*l2**2 -\
                    p2*l2*l3**3 + p6*l2*l3 + 2.0*l2 - p2*l3**4 - p10*l3**3 + p9*l3**2 + \
                        2.0*l3 - 0.1*sin_3t
qd_2 = p2*p13*l1**3*l2 - p2*p13*l1**3*l3 + p10*p13*l1**2*l2 - p10*p13*l1**2*l3  -\
        p8*p13*l1*l2**2*l3 - p10*p13*l1*l2**2 + p8*p13*l1*l2*l3**2 - p9*p13*l1*l2 + \
            p10*p13*l1*l3**2 + p9*p13*l1*l3 + p2*p13*l2**4 - p2*p13*l2**3*l3 + \
                p10*p13*l2**3  - p3*p13*l2**2*l3 - p9*p13*l2**2 + p2*p13*l2*l3**3 + \
                    p3*p13*l2*l3**2 - 2.0*p13*l2 - p2*p13*l3**4 - p10*p13*l3**3 + \
                        p9*p13*l3**2 + 2.0*p13*l3 - 0.1*cos_3t
qd_3 = -p7*l1**3 - p11*l1**2 + p10*l1*l2*l3 + p11*l1*l2 + p11*l1*l3 + l1/3 - \
    p7*l2**3 - p11*l2**2 + p11*l2*l3 + l2/3 - p7*l3**3 - p11*l3**2 + l3/3 + p12


sol = sy.solve(
    [qd_1, qd_2, qd_3],
    [l1, l2, l3]
)


print(sol)