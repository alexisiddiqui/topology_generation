import os
import sys
import subprocess
import shutil
import MDAnalysis as mda
from .pdb_modifier import run_propka, parse_pka, modify_pdb_for_ph

def get_gmx_executable():
    """Returns 'gmx_mpi' if available, otherwise 'gmx'."""
    if shutil.which("gmx_mpi"):
        return "gmx_mpi"
    return "gmx"

def run_gmx_command(arguments: list, stdin_input: str = None, cwd: str = None):
    """Runs a GROMACS command using subprocess."""
    gmx_executable = get_gmx_executable()
    command = [gmx_executable] + arguments
    
    print(f"Running command: {' '.join(command)}")
    
    # Set up environment with GMXLIB pointing to current directory
    # This ensures GROMACS can find the custom force field
    env = os.environ.copy()
    ff_base_dir = os.getcwd()
    
    # If GMXLIB is already set, append our directory; otherwise set it
    if 'GMXLIB' in env:
        env['GMXLIB'] = f"{ff_base_dir}:{env['GMXLIB']}"
    else:
        env['GMXLIB'] = ff_base_dir
    
    # Prepare stdin input - handle both string and bytes
    stdin_bytes = None
    if stdin_input is not None:
        if isinstance(stdin_input, bytes):
            stdin_bytes = stdin_input
        else:
            stdin_bytes = stdin_input.encode()
    
    result = subprocess.run(
        command,
        input=stdin_bytes,
        capture_output=True,
        text=False,  # Use binary mode for both input and output
        cwd=cwd,
        env=env,
    )
    
    if result.returncode != 0:
        print(f"Error running GROMACS command: {' '.join(command)}")
        print(f"Stdout: {result.stdout.decode()}")
        print(f"Stderr: {result.stderr.decode()}")
        raise RuntimeError("GROMACS command failed.")
    
    print(result.stdout.decode())
    print(result.stderr.decode())

def generate_topology_for_ph(pdb_file: str, ph: float, output_dir: str):
    """
    Generates a GROMACS topology for a given PDB file at a specific pH.
    """
    print(f"--- Generating topology for {pdb_file} at pH {ph} ---")

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Run PROPKA and modify PDB for target pH
    pka_file = run_propka(pdb_file)
    pka_data = parse_pka(pka_file)
    ph_pdb_file = modify_pdb_for_ph(pdb_file, pka_data, ph)
    
    base_name = os.path.basename(ph_pdb_file).replace(".pdb", "")
    
    # --- 1. pdb2gmx ---
    print("Running pdb2gmx...")
    pdb2gmx_output_gro = os.path.join(output_dir, f"{base_name}.gro")
    pdb2gmx_output_top = os.path.join(output_dir, f"{base_name}.top")
    pdb2gmx_output_itp = os.path.join(output_dir, f"{base_name}_posre.itp")

    run_gmx_command(
        ['pdb2gmx', '-f', ph_pdb_file, '-o', pdb2gmx_output_gro, '-p', pdb2gmx_output_top, '-i', pdb2gmx_output_itp, '-ff', 'amber14sb_OL21', '-water', 'tip3p', '-ignh']
    )
    print("pdb2gmx finished.")

    # --- 2. Create C-alpha position restraints ---
    print("Creating C-alpha position restraints...")
    u = mda.Universe(pdb2gmx_output_gro)
    ca_atoms = u.select_atoms("name CA")
    ca_indices = ca_atoms.indices + 1

    ca_posre_itp = os.path.join(output_dir, f"{base_name}_posre_ca.itp")
    with open(ca_posre_itp, "w") as f:
        f.write("[ position_restraints ]\n")
        f.write(";  i funct       fcx        fcy        fcz\n")
        for idx in ca_indices:
            f.write(f"{idx:4d}   1       1000       1000       1000\n")

    # Include CA restraints in the topology file
    with open(pdb2gmx_output_top, "a") as f:
        f.write(f'\n; Include C-alpha position restraints\n')
        f.write(f'#ifdef POSRES_CA\n')
        f.write(f'#include "{os.path.basename(ca_posre_itp)}"\n')
        f.write(f'#endif\n')
    print("C-alpha position restraints created.")

    # --- 3. editconf ---
    print("Running editconf...")
    editconf_output_gro = os.path.join(output_dir, f"{base_name}_box.gro")
    run_gmx_command(
        ['editconf', '-f', pdb2gmx_output_gro, '-o', editconf_output_gro, '-bt', 'cubic', '-d', '1.0']
    )
    print("editconf finished.")

    # --- 4. solvate ---
    print("Running solvate...")
    solvate_output_gro = os.path.join(output_dir, f"{base_name}_solv.gro")
    run_gmx_command(
        ['solvate', '-cp', editconf_output_gro, '-cs', 'spc216.gro', '-p', pdb2gmx_output_top, '-o', solvate_output_gro]
    )
    print("solvate finished.")

    # --- 5. Add ions ---
    print("Running grompp for ions...")
    ions_tpr = os.path.join(output_dir, "ions.tpr")
    run_gmx_command(
        ['grompp', '-f', os.path.abspath('config/ions.mdp'), '-c', solvate_output_gro, '-p', pdb2gmx_output_top, '-o', ions_tpr, '-maxwarn', '2']
    )

    print("Running genion...")
    genion_output_gro = os.path.join(output_dir, f"{base_name}_solv_ions.gro")
    run_gmx_command(
        ['genion', '-s', ions_tpr, '-p', pdb2gmx_output_top, '-o', genion_output_gro, '-neutral'],
        stdin_input="SOL\n"
    )
    print("genion finished.")

    # --- 6. Energy Minimization ---
    print("Running grompp for energy minimization...")
    em_tpr = os.path.join(output_dir, "em.tpr")
    run_gmx_command(
        ['grompp', '-f', os.path.abspath('config/minim.mdp'), '-c', genion_output_gro, '-r', pdb2gmx_output_gro, '-p', pdb2gmx_output_top, '-o', em_tpr, '-maxwarn', '2']
    )

    print("Running mdrun for energy minimization...")
    run_gmx_command(
        ['mdrun', '-v', '-deffnm', 'em'],
        cwd=output_dir
    )
    print("Energy minimization finished.")
    
    print(f"--- Topology generation for pH {ph} complete. Files are in {output_dir} ---")


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python -m scripts.protonation_pipeline.generate_topology <pdb_id>")
        sys.exit(1)
    
    pdb_id = sys.argv[1]
    
    # Assuming download_pdb is in the same package
    from .download_pdb import download_pdb
    pdb_file = download_pdb(pdb_id)

    if pdb_file:
        generate_topology_for_ph(pdb_file, ph=7.0, output_dir=f"top/{pdb_id}_pH7")
        generate_topology_for_ph(pdb_file, ph=4.0, output_dir=f"top/{pdb_id}_pH4")