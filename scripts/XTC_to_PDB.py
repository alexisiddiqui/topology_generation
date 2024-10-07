import os 

import MDAnalysis as mda 
import numpy as np


def dump_XTC(pdb_path:str, traj_path:str, frames:np.array=None, out_prefix:str=None):

    dir = os.path.dirname(traj_path)
    if out_prefix is None:
        name = os.path.basename(traj_path).split(".")[0]
    else:
        name = out_prefix

    u = mda.Universe(pdb_path, traj_path)
    print(len(u.trajectory))
    if frames is None:
        frames = np.arange(0, len(u.trajectory))

    dir = os.path.join(dir, name)
    os.system(f"rm -rf {dir}")
    os.makedirs(dir, exist_ok=True)
    
    for idx, ts in enumerate(u.trajectory[frames]):
        new_pdb_name = f"{name}_c{idx}.pdb"
        new_pdb_path = os.path.join(dir, new_pdb_name)
        u.atoms.write(new_pdb_path)


if __name__ == '__main__':

    # traj_path = "/home/alexi/Documents/topology_generation/RW_10/BPTI/BPTI_shaw_small_recl_csize2_10_.xtc"
    # pdb_path = "/home/alexi/Documents/topology_generation/RW_10/BPTI/P00974_60_1_af_sample_127_10000_protonated.pdb"
    # dump_XTC(pdb_path,traj_path, out_prefix="BPTI_10")

    traj_path = "/home/alexi/Documents/topology_generation/RW_10/BRD4/BRD4_af_small_recl_csize2_6_.xtc"
    pdb_path = "/home/alexi/Documents/topology_generation/RW_10/BRD4/BRD4_APO_484_1_af_sample_127_10000_protonated.pdb"
    dump_XTC(pdb_path,traj_path, out_prefix="BRD4_6")

    traj_path = "/home/alexi/Documents/topology_generation/RW_10/BRD4a/BRD4a_af_small_recl_csize2_6_.xtc"
    dump_XTC(pdb_path,traj_path, out_prefix="BRD4a_6")

    traj_path = "/home/alexi/Documents/topology_generation/RW_10/BRD4b/BRD4b_af_small_recl_csize2_6_.xtc"
    dump_XTC(pdb_path,traj_path, out_prefix="BRD4b_6")

    # traj_path = "/home/alexi/Documents/topology_generation/RW_10/LXR/LXRa_af_small_recl_csize2_10_.xtc"
    # pdb_path = "/home/alexi/Documents/topology_generation/RW_10/LXR/LXRa200_1_af_sample_127_10000_protonated.pdb"
    # dump_XTC(pdb_path,traj_path, out_prefix="LXR_10")

    # traj_path = "/home/alexi/Documents/topology_generation/RW_10/HOIP/HOIP_af_small_recl_csize2_10_.xtc"
    # pdb_path = "/home/alexi/Documents/topology_generation/RW_10/HOIP/HOIP_apo697_1_af_sample_127_10000_protonated.pdb"
    # dump_XTC(pdb_path,traj_path, out_prefix="HOIP_10")

    # traj_path = "/home/alexi/Documents/topology_generation/RW_10/MBP/MBP_af_small_recl_csize2_10_.xtc"
    # pdb_path = "/home/alexi/Documents/topology_generation/RW_10/MBP/MBP_wt_1_af_sample_127_10000_protonated.pdb"
    # dump_XTC(pdb_path,traj_path, out_prefix="MBP_10")




