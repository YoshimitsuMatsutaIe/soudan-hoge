"""実行例"""

import numpy as np

from simulator import Simulator


def ex_default():
    """デフォルトの実行例"""
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.reproduce_state()
    hoge.plot_basic()
    hoge.make_aniation()



def en_tuiju_by_rmp():
    """円軌道追従
    """
    
    N = 5
    TIME_SPAN = 6
    
    env_param = {
        "goal" : {
            4 : {
                "name" : "Circle",
                "param" :{
                    "r" : 0.2,
                    'center' : [0.0, 0.0, 0.75],
                    "omega" : 1.5,
                    "alpha" : 0,
                    "beta" : 0,
                    "gamma" : 0
                }
            }
        },
        "target_param" : None,
    }
    
    controller_param = {
        "name" : "rmp",
        "attractor" : {
            "max_speed" : 1800,
            "gain" : 15000,
            "a_damp_r" : 0.05,
            "sigma_W" : 1,
            "sigma_H" : 1,
            "A_damp_r" : 0.01,
        },
        "jlavoidance" : {
            "gamma_p" : 4,
            "gamma_d" : 1,
            "lam" : 1
        }
    }
    
    hoge = Simulator(N=N, controller_param=controller_param)
    
    hoge.run(TIME_SPAN)
    hoge.reproduce_state()
    hoge.save_data()
    hoge.plot_basic()
    hoge.make_aniation()



def en_tuiju_by_pdfb():
    """pdfbで円追従"""
    N = 5
    TIME_SPAN = 6
    
    env_param = {
        "goal" : {
            4 : {
                "name" : "Circle",
                "param" :{
                    "r" : 0.2,
                    'center' : [0.0, 0.0, 0.75],
                    "omega" : 1.5,
                    "alpha" : 0,
                    "beta" : 0,
                    "gamma" : 0
                }
            }
        },
        "target_param" : None,
    }
    
    controller_param = {
        "name" : "pdfb",
        "Kp" : 250,
        "Kd" : 10,
    }
    
    hoge = Simulator(N=N, controller_param=controller_param)
    
    hoge.run(TIME_SPAN)
    hoge.reproduce_state()
    hoge.save_data()
    hoge.plot_basic()
    hoge.make_aniation()


def point_3_tuiju():
    """3点を追従
    
    
    """
    N = 5
    TIME_SPAN = 6

    env_param = {
        "goal" : {
            4 : {
                "name" : "Point",
                "param" :{'center' : [0.14, 0, 0.63]}
            },
            3 : {
                "name" : "Point",
                "param" :{'center' : [0.2, 0.05, 0.5]}
            },
            2 : {
                "name" : "Point",
                "param" :{'center' : [0.14, 0, 0.4]}
            }
        },
        "target_param" : None,
    }
    
    controller_param = {
        "name" : "rmp",
        "attractor" : {
            "max_speed" : 1800,
            "gain" : 15000,
            "a_damp_r" : 0.05,
            "sigma_W" : 1,
            "sigma_H" : 1,
            "A_damp_r" : 0.01,
        },
        "jlavoidance" : {
            "gamma_p" : 4,
            "gamma_d" : 1,
            "lam" : 1
        }
    }
    
    hoge = Simulator(
        N=N,
        env_param=env_param,
        controller_param=controller_param
    )
    
    hoge.run(TIME_SPAN)
    hoge.reproduce_state()
    hoge.save_data()
    hoge.plot_basic()
    hoge.make_aniation()


# 就活
def ex_yasukawa():
    """安川電機プレゼン用
    
    ・研究概要 動画不可
    """
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.reproduce_state()
    hoge.plot_basic()
    hoge.make_aniation()


def ex_denso():
    """デンソー用
    
    ・スライド．動画可?
    """
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.reproduce_state()
    hoge.plot_basic()
    hoge.make_aniation()


def ex_kawasaki():
    """川崎重工プレゼン用
    
    
    """
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.reproduce_state()
    hoge.plot_basic()
    hoge.make_aniation()


def ex_fujikoshi():
    """不二越 研究プレゼン
    
    動画不可？
    """
    hoge = Simulator(N=5)
    
    hoge.run()
    hoge.reproduce_state()
    hoge.plot_basic()
    hoge.make_aniation()



if __name__ == '__main__':
    #en_tuiju_by_rmp()
    
    
    #en_tuiju_by_pdfb()
    
    
    
    point_3_tuiju()