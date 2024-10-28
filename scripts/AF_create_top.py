import os
import sys
import subprocess
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

def AF10K_top_gen(pdb_path:str):

    pdb_name = os.path.basename(pdb_path)


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

    subprocess.run(box_command, cwd=new_top_dir, check=True)

    # 4. Solvate the system 

    solvate_command = [gmx,
                        "solvate",
                        "-cp",
                        os.path.join(name + "_box.gro"),
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
                    os.path.join(name + ".gro"),
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
                    name + "_em"]

    subprocess.run(em_command, cwd=new_top_dir, check=True)

    # Copy the correct only the required files to the final directory

    clean_top_dir = os.path.join("clean_top", name)
    os.system("rm -rf " + clean_top_dir)
    os.makedirs(clean_top_dir, exist_ok=True)

    os.system(f"cp {new_top_dir}/{name}.top {clean_top_dir}")
    os.system(f"cp {new_top_dir}/{name}_solv_ions.gro {clean_top_dir}")
    os.system(f"cp {new_top_dir}/{name}.itp {clean_top_dir}")
    os.system(f"cp {new_top_dir}/{name}_ca.itp {clean_top_dir}")
    os.system(f"cp {new_top_dir}/{name}_em.gro {clean_top_dir}")


def iterate_over_dir(dir_path:str):
    for pdb in os.listdir(dir_path):
        if pdb.endswith(".pdb"):
            AF10K_top_gen(os.path.join(dir_path, pdb))



if __name__ == "__main__":

    if len(sys.argv) > 1:

        if sys.argv[1] == "-h":
            print("Usage: python AF_create_top.py <dir_path>")
            sys.exit(0)
        
        if isinstance(sys.argv[1], str):
            dir_path = sys.argv[1]
            iterate_over_dir(dir_path)

    else:
        protein_names = ["BPTI", "BRD4", "HOIP", "LXR", "MBP"]

        for protein in protein_names[:1]:

            base_dir_path = "/home/alexi/Documents/topology_generation/max_pLDDT_cluster10"

            path_suffix = f"{protein}_10"

            dir_path = os.path.join(base_dir_path, protein, path_suffix)

            iterate_over_dir(dir_path)