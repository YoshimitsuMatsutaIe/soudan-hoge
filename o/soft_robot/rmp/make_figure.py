"""グラフを作る"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "DejaVu Serif"





# def load_data_from_csv(data_path):
#     actu_data = np.loadtxt(data_path + "/actuator.csv", delimiter=',', dtype='float64')
#     task_data = np.loadtxt(data_path + "/task.csv", delimiter=',', dtype='float64')
    
#     _, n = actu_data.shape
#     N = (n-1) // 6
    
#     _, n = task_data.shape
#     T = (n-1-N*3)
    
#     return actu_data, task_data, N, T


def load_df_from_csv(data_path):
    actu_df = pd.read_csv(data_path + "/actuator.csv")
    task_df = pd.read_csv(data_path + "/task.csv")
    
    n = len(actu_df.columns)
    N = (n-1) // 6
    
    n = len(task_df.columns)
    T = (n-1-N*3) // 4
    
    return actu_df, task_df, N, T



# def plot_error(pd, cd=None, N, T):
#     """誤差をプロット"""
    
#     figs = []
#     for i in range(T):
#         figs.append(plt.figure())
#         ax = figs[i].add_subplot(111)
#         ax.plot(
#             pd[:, 0], pd[:, 1+3*N+i],
#             label="RMP", linestyle = 'solid'
#         )
        
#         if cd is not None:
#             ax.plot(
#                 cd[:, 0], cd[:, 1+3*N+i],
#                 label="PD-fb", linestyle = 'dashed'
#             )
        
#         ax.set_xlabel(r'Time $\it{t}$ [s]')
#         ax.set_ylabel(r"Position Error [m]")
#         ax.set_xlim(0, pd[-1,0])
        
#         ax.legend()
#         ax.grid()
        
#         figs[i].savefig('error_' + )
    
#     plt.show()



if __name__ == "__main__":
    
    # 比較するデータフォルダの絶対or相対パス
    pdfb_path = r"/home/matsuta/src/ctrlab2021_soudan/result/2022-01-28--19-10-50"
    rmp_path = r"/home/matsuta/src/ctrlab2021_soudan/result/2022-01-28--19-10-50"

    pad, ptd, N, T = load_df_from_csv(rmp_path)
    print(pad)