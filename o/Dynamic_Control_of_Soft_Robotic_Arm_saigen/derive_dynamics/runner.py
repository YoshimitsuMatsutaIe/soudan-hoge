import sympy as sy

import dynamics


if __name__ == ("__main__"):
    N = 3
    hoge = dynamics.Dynamics(3)
    
    
    for i, Phi in enumerate(hoge.Phi_s):
        name = "Phi" + str(i) + ".txt"
        f = open(name, 'w')
        f.write(str(Phi))
        f.close()
    
    for i, Theta in enumerate(hoge.Theta_s):
        name = "Theta" + str(i) + ".txt"
        f = open(name, 'w')
        f.write(str(Theta))
        f.close()
    
    
    f = open("M.txt", 'w')
    f.write(str(hoge.M))
    f.close()

    f = open("C.txt", 'w')
    f.write(str(hoge.C))
    f.close()
    
    f = open("G.txt", 'w')
    f.write(str(hoge.G))
    f.close()
