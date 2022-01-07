import sympy as sy
from sympy.printing.numpy import NumPyPrinter
from sympy import julia_code

from sympy.utilities.codegen import codegen

import dynamics
import utils
import os



def make_function():
    """シミュレーションに必要な式を生成"""
    
    N = 1  # セクションの数
    
    
    cwd = os.path.dirname(__file__)
    base = cwd + "/../derived/N_is_" + str(N)
    
    
    # 時間がかかるので細かくセーブしながら実行
    # どこまで終わってるか確認しながら手作業で実行
    
    dir_name = base + '/obj'
    os.makedirs(dir_name, exist_ok=True)
    
    
    if not os.path.isfile(dir_name + "/until_kinematics" + ".binaryfile"):
        hoge = dynamics.Dynamics(N)
        utils.save_obj_by_picke(hoge, dir_name, "/until_kinematics",)
    
    
    if not os.path.isfile(dir_name + "/until_M_omega" + ".binaryfile"):
        hoge = utils.load_obj_from_picle(dir_name, "/until_kinematics",)
        hoge.set_M_omega()
        utils.save_obj_by_picke(hoge, dir_name, "/until_M_omega")
    
    
    if not os.path.isfile(dir_name + "/until_M_omega_dot" + ".binaryfile"):
        hoge = utils.load_obj_from_picle(dir_name, "/until_M_omega",)
        hoge.set_M_omega_dot()
        utils.save_obj_by_picke(hoge, dir_name, "/until_M_omega_dot")
    
    
    if not os.path.isfile(dir_name + "/until_M_v" + ".binaryfile"):
        hoge = utils.load_obj_from_picle(dir_name, "/until_M_omega_dot",)
        hoge.set_M_v()
        utils.save_obj_by_picke(hoge, dir_name, "/until_M_v")
    
    
    if not os.path.isfile(dir_name + "/until_M_v_dot" + ".binaryfile"):
        hoge = utils.load_obj_from_picle(dir_name, "/until_M_v",)
        hoge.set_M_v_dot()
        utils.save_obj_by_picke(hoge, dir_name, "/until_M_v_dot")
    
    
    if not os.path.isfile(dir_name + "/until_M" + ".binaryfile"):
        hoge = utils.load_obj_from_picle(dir_name, "/until_M_v_dot",)
        hoge.set_M()
        utils.save_obj_by_picke(hoge, dir_name, "/until_M")
    
    
    if not os.path.isfile(dir_name + "/until_M_dot" + ".binaryfile"):
        hoge = utils.load_obj_from_picle(dir_name, "/until_M",)
        hoge.set_M_dot()
        utils.save_obj_by_picke(hoge, dir_name, "/until_M_dot")
    
    
    if not os.path.isfile(dir_name + "/until_C" + ".binaryfile"):
        hoge = utils.load_obj_from_picle(dir_name, "/until_M_dot",)
        hoge.set_C()
        utils.save_obj_by_picke(hoge, dir_name, "/until_C")
    
    
    if not os.path.isfile(dir_name + "/until_G_g" + ".binaryfile"):
        hoge = utils.load_obj_from_picle(dir_name, "/until_C",)
        hoge.set_G_g()
        utils.save_obj_by_picke(hoge, dir_name, "/until_G_g")
    
    
    if not os.path.isfile(dir_name + "/until_G_e" + ".binaryfile"):
        hoge = utils.load_obj_from_picle(dir_name, "/until_G_g",)
        hoge.set_G_e()
        utils.save_obj_by_picke(hoge, dir_name, "/until_G_e")
    
    
    if not os.path.isfile(dir_name + "/until_G" + ".binaryfile"):
        hoge = utils.load_obj_from_picle(dir_name, "/until_G_e",)
        hoge.set_G()
        utils.save_obj_by_picke(hoge, dir_name, "/until_G")
    else:
        hoge = utils.load_obj_from_picle(dir_name, "/until_G",)


    # # 数式を整理（やばいくらい時間かかる）
    # print("数式整理中...")
    # for i, Phi in enumerate(hoge.Phi_s):
    #     hoge.Phi_s[i] = sy.simplify(Phi)
    # for i, Theta in enumerate(hoge.Theta_s):
    #     hoge.Theta_s[i] = sy.simplify(Theta)
    # hoge.M = sy.simplify(hoge.M)
    # hoge.C = sy.simplify(hoge.C)
    # hoge.G = sy.simplify(hoge.G)

    # 数式を生成
    dir_name = base + '/eqs/txt_style'
    os.makedirs(dir_name, exist_ok=True)
    
    
    for i, Phi in enumerate(hoge.Phi_s):
        name = dir_name + "/Phi" + str(i) + ".txt"
        f = open(name, 'w')
        f.write(str(Phi))
        f.close()
    
    for i, Theta in enumerate(hoge.Theta_s):
        name = dir_name + "/Theta" + str(i) + ".txt"
        f = open(name, 'w')
        f.write(str(Theta))
        f.close()
    
    for i in range(3*N):
        for j in range(3*N):
            f = open(dir_name + "/M" + str(i) + "_" + str(j) + ".txt", 'w')
            f.write(str(hoge.M[i, j]))
            f.close()

    for i in range(3*N):
        for j in range(3*N):
            f = open(dir_name + "/C" + str(i) + "_" + str(j) + ".txt", 'w')
            f.write(str(hoge.C[i, j]))
            f.close()
    
    for i in range(3*N):
        f = open(dir_name + "/G" + str(i) + ".txt", 'w')
        f.write(str(hoge.G[i, :]))
        f.close()

    
    # numpyスタイルの数式を生成
    dir_name = base + '/eqs/numpy_style'
    os.makedirs(dir_name, exist_ok=True)
    
    numpy_word = "import numpy\ndef f(q_large, xi_large, q_dot_large):\n    return "
    
    for i, Phi in enumerate(hoge.Phi_s):
        name = dir_name + "/Phi" + str(i) + ".py"
        f = open(name, 'w')
        f.write(numpy_word)
        f.write(NumPyPrinter().doprint(Phi))
        f.close()
    
    for i, Theta in enumerate(hoge.Theta_s):
        name = dir_name + "/Theta" + str(i) + ".py"
        f = open(name, 'w')
        f.write(numpy_word)
        f.write(NumPyPrinter().doprint(Theta))
        f.close()
    
    for i in range(3*N):
        for j in range(3*N):
            f = open(dir_name + "/M" + str(i) + "_" + str(j) + ".py", 'w')
            f.write(numpy_word)
            f.write(NumPyPrinter().doprint(hoge.M[i, j]))
            f.close()

    for i in range(3*N):
        for j in range(3*N):
            f = open(dir_name + "/C" + str(i) + "_" + str(j) + ".py", 'w')
            f.write(numpy_word)
            f.write(NumPyPrinter().doprint(hoge.C[i, j]))
            f.close()
    
    for i in range(3*N):
        f = open(dir_name + "/G" + str(i) + ".py", 'w')
        f.write(numpy_word)
        f.write(NumPyPrinter().doprint(hoge.G[i, 0]))
        f.close()
    

    # じゅりあのコード生成
    dir_name = base + '/eqs/julia_style'
    os.makedirs(dir_name, exist_ok=True)
    
    julia_word = "function f(q_large::Vector{T}, xi_large::Vector{T}, q_dot_large::Vector{T}) where T\n    "
    
    for i, Phi in enumerate(hoge.Phi_s):
        name = dir_name + "/Phi" + str(i) + ".jl"
        f = open(name, 'w')
        f.write("module " + "Phi" + str(i) + "\n")
        f.write(julia_word)
        f.write(julia_code(Phi))
        f.write("\nend\nend")
        f.close()
    
    for i, Theta in enumerate(hoge.Theta_s):
        name = dir_name + "/Theta" + str(i) + ".jl"
        f = open(name, 'w')
        f.write("module " + "Theta" + str(i) + "\n")
        f.write(julia_word)
        f.write(julia_code(Theta))
        f.write("\nend\nend")
        f.close()
    
    for i in range(3*N):
        for j in range(3*N):
            f = open(dir_name + "/M" + str(i) + "_" + str(j) + ".jl", 'w')
            f.write("module " + "M" + str(i) + "_" + str(j) + "\n")
            f.write(julia_word)
            f.write(julia_code(hoge.M[i, j]))
            f.write("\nend\nend")
            f.close()

    for i in range(3*N):
        for j in range(3*N):
            f = open(dir_name + "/C" + str(i) + "_" + str(j) + ".jl", 'w')
            f.write("module " + "C" + str(i) + "_" + str(j) + "\n")
            f.write(julia_word)
            f.write(julia_code(hoge.C[i, j]))
            f.write("\nend\nend")
            f.close()
    
    for i in range(3*N):
        f = open(dir_name + "/G" + str(i) + ".jl", 'w')
        f.write("module " + "G" + str(i) + "\n")
        f.write(julia_word)
        f.write(julia_code(hoge.G[i, 0]))
        f.write("\nend\nend")
        f.close()
    


    # おまけでCのコード生成
    

    
    dir_name = base + '/eqs/c_src/'
    os.makedirs(dir_name, exist_ok=True)
    
    
    def gen_c(f, name, dir_name):
        [(c_name, c_code), (h_name, c_header)] = codegen(
            name_expr=(name, f),
            language="C",
            project= name + "project",
            to_files=False
        )
        
        f = open(dir_name + c_name, 'w')
        f.write(c_code)
        f.close()

        f = open(dir_name + h_name, 'w')
        f.write(c_header)
        f.close()
    
    
    for i, Phi in enumerate(hoge.Phi_s):
        name = "Phi" + str(i)
        gen_c(Phi, name, dir_name)
    
    for i, Theta in enumerate(hoge.Theta_s):
        name = "Theta" + str(i)
        gen_c(Theta, name, dir_name)
    
    
    gen_c(hoge.M, "M", dir_name)
    gen_c(hoge.C, "C", dir_name)
    gen_c(hoge.G, "G", dir_name)


if __name__ == ("__main__"):
    make_function()