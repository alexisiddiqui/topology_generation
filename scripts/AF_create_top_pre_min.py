import os
import sys
import subprocess
from pytrr import GroTrrReader
import pyedr

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
# os.chdir("/home/alexi/Documents/topology_generation/")
# # downlaod the AMBER FF14SB force fid
# amber_14_url = "https://ftp.gromacs.org/contrib/forcefields/amber14sb_OL15.ff_corrected-Na-cation-params.tar.gz"
# amber_14_url = "https://fch.upol.cz/ff_ol/amber14sb_OL21.ff.tar.gz"

# subprocess.run(["wget", amber_14_url])
# # subprocess.run(["tar", "-xvf", "amber14sb_OL15.ff_corrected-Na-cation-params.tar.gz", "-C", os.getcwd()])
# subprocess.run(["tar", "-xvf", "amber14sb_OL21.ff.tar.gz", "-C", os.getcwd()])
# subprocess.run(["rm  amber14sb*.tar.gz"], shell=True)

gmx = "gmx_mpi"
os.environ["GMXLIB"] = os.path.join(os.getcwd())




# This script is used to generate topology files from PDB files from localcolabfold sampling

# point to the directory where the pdb files are stored - generate topology files for each PDB files in the directory
# Follow standard GROMACs topology generation procedure with AMBER FF14SB force field and TIP3P water model
# Main modification to the procedure is to create a restraints file for the protein CA atoms using the original structure during equilabration 

# 1. Generate topology files using pdb2gmx

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

#     assert len(vacuo_trr) == len(solvent_trr) == len(names), "The number of entries in the two trr arrays must be the same"

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
        # print(vacuo_trr[idx].astype(float))
        # print(solvent_trr[idx].astype(float))
    #     if not isinstance(vacuo_trr[idx], np.ndarray):
    #         vacuo_trr[idx] = np.array(vacuo_trr[idx])
    #     if not isinstance(solvent_trr[idx], np.ndarray):
    #         solvent_trr[idx] = np.array(solvent_trr[idx]))
        ax1.plot(vacuo_trr[idx], label=f"{name} {idx} vacuo")
        ax2.plot(solvent_trr[idx], label=f"{name} {idx} solvent")

    ax1.set_title("Vacuo")
    ax1.set_xlabel("Frame")
    ax1.set_ylabel(ylabel)
    ax1.legend()

    ax2.set_title("Solvent")
    ax2.set_xlabel("Frame")
    ax2.set_ylabel(ylabel)
    ax2.legend()

    plt.tight_layout()

    if save_path is None:
        save_path = os.getcwd()

    try:
        _names = []
        for idx, name in enumerate(names):
            name = name.split("_")[0] + str(idx)
            _names.append(name)
        name_str = "_".join(_names)
        save_path = os.path.join(save_path, f"{name_str}_min_{ylabel}_vac-sol.png")
        plt.savefig(save_path, dpi=300)

    except:
        names_str = _names[0] + "_n" + str(len(_names))
        save_path = os.path.join(save_path, f"{names_str}_min_{ylabel}_vac-sol.png")
        plt.savefig(save_path, dpi=300)

    plt.close(fig)

def AF10K_top_gen(pdb_path:str,box_size):

    pdb_name = os.path.basename(pdb_path)

    dir_path = os.path.dirname(pdb_path)

    name = pdb_name.split(".")[0]
        
    new_top_dir = os.path.join("top", name)
    os.system("rm -rf " + new_top_dir)

    os.makedirs(new_top_dir, exist_ok=True)

    dummy_posres = "121212"
    
    amber_14_dir = "amber14sb_OL21" # amber14sb_OL15

    pdb2gmx_command = [gmx,
                    "pdb2gmx",
                        "-f",
                        os.path.join(dir_path, pdb_name),
                        "-o",
                        os.path.join(name + ".gro"),
                        "-p",
                        os.path.join(name + ".top"),
                        "-i",
                        os.path.join(name + ".itp"),
                        "-water",
                        "tip3p",
                        "-ff",
                        amber_14_dir,
                        "-ignh",
                        "-posrefc",
                        dummy_posres]

    print(" ".join(pdb2gmx_command))
    #set GMXLIB

    subprocess.run(pdb2gmx_command, cwd=new_top_dir, check=True)
    # 2. Create an itp file for the protein CA atoms - use the gro file from the previous step as the reference structure

    itp_path = os.path.join(new_top_dir, name + ".itp")
    gro_path = os.path.join(new_top_dir, name + ".gro")
    import MDAnalysis as mda

    u = mda.Universe(gro_path)
    print(u.atoms.indices)
    CA_atoms = u.select_atoms("name CA")

    CA_indices = CA_atoms.indices+1
    print(CA_indices)



    heavy_POSRES_flag = "POSRES_FC"

    with open(itp_path, "r") as f:
        lines = f.readlines()
        # replace the dummy posres with the custom posres
        lines = [line.replace(dummy_posres, heavy_POSRES_flag) for line in lines]
    
    with open(itp_path, "w") as f:
        f.writelines(lines)

    custom_posres_flag = "POSRESca_FC"

    # with open(itp_path, "r") as f:
    #     lines = f.readlines()

    new_lines = []  
    for line in lines:
        if len(line.split()) > 1:
            print(line)
            atom_num = (line.split()[0])
            if atom_num.isdigit() and int(atom_num) not in CA_indices:
                continue
            else:
                # replace dummy posres with custom posres
                line = line.replace(heavy_POSRES_flag, custom_posres_flag)
                print(line)
        new_lines.append(line)

    new_itp_path = os.path.join(new_top_dir, name + "_ca.itp")
    with open(new_itp_path, "w") as f:
        f.writelines(new_lines)


    # add to the topology file

    top_path = os.path.join(new_top_dir, name + ".top")

    CA_POSRES_lines = f"""

#ifdef POSRESca
#include "{name}_ca.itp"
#endif

"""

    POSRES_flag = "#ifdef POSRES"

    with open(top_path, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if POSRES_flag in line:
            lines.insert(i+3, CA_POSRES_lines)
            break

    with open(top_path, "w") as f:
        f.writelines(lines)

    # 3. Create Box
    if box_size is None:

        box_command = [gmx,
                    "editconf",
                    "-f",
                    os.path.join(name + ".gro"),
                    "-o",
                    os.path.join(name + "_box.gro"),
                    "-bt",
                    "cubic",
                    "-d",
                    "1.0"]

    else:
        box_command = [gmx,
                    "editconf",
                    "-f",
                    os.path.join(name + ".gro"),
                    "-o",
                    os.path.join(name + "_box.gro"),
                    "-bt",
                    "cubic",
                    "-box",
                    f"{box_size}",
                    f"{box_size}",
                    f"{box_size}"]

    subprocess.run(box_command, cwd=new_top_dir, check=True)

    em_command = [gmx,
                    "grompp",
                    "-f",
                    "/home/alexi/Documents/topology_generation/config/1_minim.mdp",
                    "-c",
                    os.path.join(name + "_box.gro"),
                    "-r",
                    os.path.join(name + "_box.gro"),
                    "-p",
                    os.path.join(name + ".top"),
                    "-maxwarn",
                    "2",
                    "-o",
                    "em.tpr"]

    subprocess.run(em_command, cwd=new_top_dir, check=True)

    em_command = [gmx,
                    "mdrun",
                    "-s",
                    "em.tpr",
                    "-v",
                    "-deffnm",
                    name + "_box_em",
                    "-ntomp", "20"]


    subprocess.run(em_command, cwd=new_top_dir, check=True)

    vac_trr_path = os.path.join(new_top_dir, name + "_box_em.trr")
    vac_trr_data = read_min_energy(vac_trr_path)
    vac_edr_path = vac_trr_path.replace(".trr",".edr")



    # 4. Solvate the system 

    solvate_command = [gmx,
                        "solvate",
                        "-cp",
                        os.path.join(name + "_box_em.gro"),
                        "-cs",
                        "spc216.gro",
                        "-o",
                        os.path.join(name + "_solv.gro"),
                        "-p",
                        os.path.join(name + ".top")]


    subprocess.run(solvate_command, cwd=new_top_dir, check=True)


    # Add ions

    genion_command = [gmx,
                        "grompp",
                        "-f",
                        "/home/alexi/Documents/topology_generation/config/ions.mdp",
                        "-c",
                        os.path.join(name + "_solv.gro"),
                        "-p",
                        os.path.join(name + ".top"),
                        "-o",
                        "ions.tpr"]

    subprocess.run(genion_command, cwd=new_top_dir, check=True)


    genion_command = [gmx,
                        "genion",
                        "-s",
                        "ions.tpr",
                        "-o",
                        os.path.join(name + "_solv_ions.gro"),
                        "-p",
                        os.path.join(name + ".top"),
                        "-pname",
                        "NA",
                        "-nname",
                        "CL",
                        "-neutral",
                        "yes"]

    subprocess.run(genion_command, cwd=new_top_dir, input=b"13\n", check=True)

    # 4. Energy minimization use original structure as reference for CA atoms

    em_command = [gmx,
                    "grompp",
                    "-f",
                    "/home/alexi/Documents/topology_generation/config/1_minim.mdp",
                    "-c",
                    os.path.join(name + "_solv_ions.gro"),
                    "-r",
                    os.path.join(name + "_solv_ions.gro"),
                    "-p",
                    os.path.join(name + ".top"),
                    "-maxwarn",
                    "2",
                    "-o",
                    "em.tpr"]

    subprocess.run(em_command, cwd=new_top_dir, check=True)

    em_command = [gmx,
                    "mdrun",
                    "-s",
                    "em.tpr",
                    "-v",
                    "-deffnm",
                    name + "_em",
                    "-ntomp", "20"]

    subprocess.run(em_command, cwd=new_top_dir, check=True)

    # Copy the correct only the required files to the final directory

    clean_top_dir = os.path.join("clean_top", name)
    os.system("rm -rf " + clean_top_dir)
    try:
        os.makedirs(clean_top_dir)
    except:
        raise Exception("Could not create clean top directory - check for identical names in the directory")

    os.system(f"cp {new_top_dir}/{name}.top {clean_top_dir}")
    os.system(f"cp {new_top_dir}/{name}_solv_ions.gro {clean_top_dir}")
    os.system(f"cp {new_top_dir}/{name}.itp {clean_top_dir}")
    os.system(f"cp {new_top_dir}/{name}_ca.itp {clean_top_dir}")
    # os.system(f"cp {new_top_dir}/{name}_box_em.gro {clean_top_dir}")
    sol_trr_path = os.path.join(new_top_dir, name + "_em.trr")
    sol_trr_data = read_min_energy(sol_trr_path)


    sol_edr_path = sol_trr_path.replace(".trr",".edr")

    edr_key = "Potential"
    vac_edr_data = read_EDR_file(vac_edr_path,edr_key)
    sol_edr_data = read_EDR_file(sol_edr_path,edr_key)


    return vac_trr_data, sol_trr_data, name, vac_edr_data, sol_edr_data


def iterate_over_dir(dir_path:str, box_size=None):
    vac_trr_data, sol_trr_data, names = [],[], []
    vac_edr_data, sol_edr_data = [],[]
    for pdb in os.listdir(dir_path):
        if pdb.endswith(".pdb"):
            vac_trr, sol_trr, name, vac_edr, sol_edr = AF10K_top_gen(os.path.join(dir_path, pdb),box_size)
            vac_trr_data.append(vac_trr)
            sol_trr_data.append(sol_trr)
            vac_edr_data.append(vac_edr)
            sol_edr_data.append(sol_edr)

            names.append(name)


    _names = [name.split("_")[0] for name in names]


    # clip really large values in trr data to 10000
    for idx, _ in enumerate(names):
        vac_trr_data[idx][vac_trr_data[idx] > 10000] = 10000.1
        sol_trr_data[idx][sol_trr_data[idx] > 10000] = 10000.1
        vac_edr_data[idx][vac_edr_data[idx] > 10000] = 10000.1
        sol_edr_data[idx][sol_edr_data[idx] > 10000] = 10000.1




    # create dataframe of trr data with vac and sol for each frame, with names
    data = {"name": names, "vac": vac_trr_data, "sol": sol_trr_data,  "e_vac": vac_edr_data, "e_sol": sol_edr_data}

    csv_name = "_".join(_names) + "_min_force_energy.csv"
    csv_path = os.path.join(dir_path, csv_name)

    df = pd.DataFrame(data)

    df.to_pickle(csv_path)  



    plot_min_energy(vac_trr_data, sol_trr_data, _names, ylabel="Force", save_path=dir_path)
    # df.to_csv(csv_path, index=False, float_format='%.2f')

    plot_min_energy(vac_edr_data, sol_edr_data, _names, save_path=dir_path)





if __name__ == "__main__":

    # set GMXLIB to include the AMBER FF14SB force field
    # os.environ["GMXLIB"] = os.path.join(os.getcwd())

    # check that file is being run in the correct directory
    assert os.path.basename(os.getcwd()) == "topology_generation", "Please run this script in the topology_generation directory"


    if len(sys.argv) > 1:

        if sys.argv[1] == "-h":
            print("Usage: python AF_create_top.py <dir_path>")
            sys.exit(0)
        
        if isinstance(sys.argv[1], str):
            dir_path = sys.argv[1]
            iterate_over_dir(dir_path)

    else:

        # protein_names = ["BPTI" ,"BRD4",  "LXR", "MBP"]
        # protein_names = ["BPTI"]

        # dir_path = "/home/alexi/Documents/topology_generation/max_pLDDT/"
        # # dir_path = "/home/alexi/Documents/topology_generation/MBP"
        # for pt in protein_names:
        #     path = os.path.join(dir_path, pt)
        #     print(path)
        #     iterate_over_dir(path)

        protein_names = ["HOIP"]

        dir_path = "/home/alexi/Documents/topology_generation/max_pLDDT/"
        for pt in protein_names:
            path = os.path.join(dir_path, pt)
            print(path)
            iterate_over_dir(path, box_size=15.0)

        