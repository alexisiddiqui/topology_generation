
import os 
import numpy as np
from pytrr import GroTrrReader
import matplotlib.pyplot as plt
import matplotlib
# set backend for vscode small plot
import pandas as pd

import seaborn as sns            
matplotlib.use('TkAgg')

import pyedr



def read_min_energy(trr_path, key:str='f'):


    with GroTrrReader(trr_path) as trrfile:

         # get length of trr

        trr_data = []

        for idx, frame in enumerate(trrfile):
            # print(frame['f_size'])
            frame_data = trrfile.get_data()
            print(frame_data[key].sum())
            trr_data.append(frame_data[key].sum())

        return np.array(trr_data)
    

def read_EDR_file(edr_path, key:str="Potential"):

    edr_data = pyedr.edr_to_dict(edr_path)

    print(edr_data.keys())

    e_pot = edr_data[key]

    return e_pot




# def plot_min_energy(vacuo_trr:list[np.ndarray], solvent_trr:list[np.ndarray], names: list[str]=None, save_path:str=None):

#     if names is None:
#         names = list(range(1,len(vacuo_trr)+1))
#         # convert to str
#     names = [str(name) for name in names]

#     # assert len(vacuo_trr) == len(solvent_trr) == len(names), "The number of entries in the two trr arrays must be the same"

#     # plot two subplots - plot the vacuo and solvent trr data
#     plt.figure(figsize=(10, 5))

#     # create two subplots
#     ax1 = plt.subplot(1, 2, 1)
#     ax2 = plt.subplot(1, 2, 2)


#     # label each entry by name
#     for idx, name in enumerate(names):
#         ax1.plot(vacuo_trr[idx], label=f"{name} vacuo")
#         ax2.plot(solvent_trr[idx], label=f"{name} solvent")


#     ax1.legend()
#     ax2.legend()

#     if save_path is None:
#         name_str = "_".join(names)
#         save_path = f"{name_str}_min_energy_vac-sol.png"



#     try:
#         save_path = f"{name_str}_min_energy_vac-sol.png"
#         plt.savefig(save_path, dpi=300)
    
#     except:
#         _names = []
#         for idx,name in enumerate(names):
#             name = name.split("_")[0]+str(idx)
#             _names.append(name)
#         name_str = "_".join(names)
#         save_path = f"{name_str}_min_energy_vac-sol.png"
#         plt.savefig(save_path, dpi=300)

#     # plt.show()
#     plt.close()

def plot_min_energy(vacuo_trr: list[np.ndarray], solvent_trr: list[np.ndarray], names: list[str] = None, save_path: str = None, ylabel="Energy"):
    if names is None:
        names = list(range(1, len(vacuo_trr) + 1))

    # convert to str
    names = [str(name) for name in names]
    # print(vacuo_trr.shape)
    # print(solvent_trr.shape)
    # assert len(vacuo_trr) == len(solvent_trr) == len(names), "The number of entries in the two trr arrays must be the same"

    # plot two subplots - plot the vacuo and solvent trr data
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    # label each entry by name
    for idx, name in enumerate(names):
        print(name)
        print(vacuo_trr[idx].astype(float))
        print(solvent_trr[idx].astype(float))
    #     if not isinstance(vacuo_trr[idx], np.ndarray):
    #         vacuo_trr[idx] = np.array(vacuo_trr[idx])
    #     if not isinstance(solvent_trr[idx], np.ndarray):
    #         solvent_trr[idx] = np.array(solvent_trr[idx]))
        ax1.plot(vacuo_trr[idx], label=f"{name} vacuo")
        ax2.plot(solvent_trr[idx], label=f"{name} solvent")

    ax1.set_title("Vacuo")
    ax1.set_xlabel("Frame")
    ax1.set_ylabel("Energy")
    ax1.legend()

    ax2.set_title("Solvent")
    ax2.set_xlabel("Frame")
    ax2.set_ylabel("Energy")
    ax2.legend()

    plt.tight_layout()

    if save_path is None:
        name_str = "_".join(names)
        save_path = f"{name_str}_min_energy_vac-sol.png"

    try:
        plt.savefig(save_path, dpi=300)
    except:
        _names = []
        for idx, name in enumerate(names):
            name = name.split("_")[0] + str(idx)
            _names.append(name)
        name_str = "_".join(_names)
        save_path = f"{name_str}_min_energy_vac-sol.png"
        plt.savefig(save_path, dpi=300)

    plt.close(fig)



if __name__ == '__main__':
    # import ast

    vac_trr_path = "/home/alexi/Documents/topology_generation/top/P00974_60_1_af_sample_127_10000_protonated/P00974_60_1_af_sample_127_10000_protonated_box_em.trr"
    # # trr_path = "/home/alexi/Documents/topology_generation/top/BRD4_APO_484_1_af_sample_127_10000_protonated/BRD4_APO_484_1_af_sample_127_10000_protonated_em.trr"
    # vac_trr_data = read_min_energy(trr_path)
    # # print(vac_trr_data)
    vac_edr_path = vac_trr_path.replace(".trr",".edr")

    sol_trr_path = "/home/alexi/Documents/topology_generation/top/P00974_60_1_af_sample_127_10000_protonated/P00974_60_1_af_sample_127_10000_protonated_em.trr"
    # # trr_path = "/home/alexi/Documents/topology_generation/top/BRD4_APO_484_1_af_sample_127_10000_protonated/BRD4_APO_484_1_af_sample_127_10000_protonated_em.trr"
    # sol_trr_data = read_min_energy(trr_path)
    # # print(sol_trr_data)
    sol_edr_path = sol_trr_path.replace(".trr",".edr")

    edr_key = "Potential"
    vac_edr_data = read_EDR_file(vac_edr_path,edr_key)
    sol_edr_data = read_EDR_file(sol_edr_path,edr_key)

    print(vac_edr_data)
    print(sol_edr_path)

    # plot_min_energy([vac_trr_data], [sol_trr_data])

    csv_path = "/home/alexi/Documents/topology_generation/test/HOIP_min_energy.pkl"
    df = pd.read_pickle(csv_path)
    # print(df)
    # vac_trr_data = df['vac'].to_num/py().astype(float)
    # print(vac_trr_data)
    # vac_trr_data = ast.literal_eval(vac_trr_data)


    # sol_trr_data = df['sol'].to_numpy(). .astype(float)
    

    names = df['name'].values
    _names = [name.split("_")[0] for name in names]

    vac_trr_data = df['vac'].values
    # convert str to float
    sol_trr_data = df['sol'].values
    # import ast

    # for idx, name in enumerate(names):
    #     vac_trr_data[idx] = np.array(ast.literal_eval(vac_trr_data[idx]))
    #     sol_trr_data[idx] = np.array(ast.literal_eval(sol_trr_data[idx]))


    plot_min_energy([vac_edr_data], [sol_edr_data], names=[_names[0]])