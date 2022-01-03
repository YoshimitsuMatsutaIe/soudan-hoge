import sympy as sy

import dynamics
import utils
import os



def make_function():
    """シミュレーションに必要な式を生成"""
    
    N = 2  # セクションの数
    
    
    cwd = os.path.dirname(__file__)
    base = cwd + "/result/N_is_" + str(N)
    
    
    # 時間がかかるので細かくセーブしながら実行
    # どこまで終わってるか確認しながら手作業で実行
    
    dir_name = base + '/obj'
    os.makedirs(dir_name, exist_ok=True)
    
    
    hoge = dynamics.Dynamics(N)
    hoge = utils.save_obj_by_picke(hoge, dir_name, "/until_kinematics",)
    
    
    hoge = utils.load_obj_from_picle(dir_name, "/until_kinematics",)
    hoge.set_M_omega()
    utils.save_obj_by_picke(hoge, dir_name, "/until_M_omega")
    
    
    hoge = utils.load_obj_from_picle(dir_name, "/until_M_omega",)
    hoge.set_M_omega_dot()
    utils.save_obj_by_picke(hoge, dir_name, "/until_M_omega_dot")
    
    
    hoge = utils.load_obj_from_picle(dir_name, "/until_M_omega_dot",)
    hoge.set_M_v()
    utils.save_obj_by_picke(hoge, dir_name, "/until_M_v")
    
    
    hoge = utils.load_obj_from_picle(dir_name, "/until_M_v",)
    hoge.set_M_v_dot()
    utils.save_obj_by_picke(hoge, dir_name, "/until_M_v_dot")
    
    
    hoge = utils.load_obj_from_picle(dir_name, "/until_M_v_dot",)
    hoge.set_M()
    utils.save_obj_by_picke(hoge, dir_name, "/until_M")
    
    
    hoge = utils.load_obj_from_picle(dir_name, "/until_M",)
    hoge.set_M_dot()
    utils.save_obj_by_picke(hoge, dir_name, "/until_M_dot")
    
    
    hoge = utils.load_obj_from_picle(dir_name, "/until_M_dot",)
    hoge.set_C()
    utils.save_obj_by_picke(hoge, dir_name, "/until_C")
    
    
    hoge = utils.load_obj_from_picle(dir_name, "/until_C",)
    hoge.set_G_g()
    utils.save_obj_by_picke(hoge, dir_name, "/until_G_g")
    
    
    hoge = utils.load_obj_from_picle(dir_name, "/until_G_g",)
    hoge.set_G_e()
    utils.save_obj_by_picke(hoge, dir_name, "/until_G_e")
    
    
    hoge = utils.load_obj_from_picle(dir_name, "/until_G_e",)
    hoge.set_G()
    utils.save_obj_by_picke(hoge, dir_name, "/until_G")
    


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
    
    
    f = open(dir_name + "/M.txt", 'w')
    f.write(str(hoge.M))
    f.close()

    f = open(dir_name + "/C.txt", 'w')
    f.write(str(hoge.C))
    f.close()
    
    f = open(dir_name + "/G.txt", 'w')
    f.write(str(hoge.G))
    f.close()




if __name__ == ("__main__"):
    make_function()