import os
import sys
import subprocess
import shutil
from typing import List, Dict

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

def run_equilibration_step(
    step_num: int,
    mdp_file: str,
    input_gro: str,
    input_top: str,
    input_cpt: str,
    output_dir: str,
    base_name: str,
    reference_gro: str = None  # Add this parameter
) -> Dict[str, str]:
    """Runs a single step of equilibration."""
    print(f"--- Running equilibration step {step_num} ---")
    
    step_name = f"step{step_num}_equil"
    output_tpr = os.path.join(output_dir, f"{step_name}.tpr")
    
    grompp_args = ['grompp', '-maxwarn', '2', '-f', mdp_file, '-c', input_gro, '-p', input_top, '-o', output_tpr]
    
    # Add reference structure for position restraints
    if reference_gro:
        grompp_args.extend(['-r', reference_gro])
    
    if input_cpt:
        grompp_args.extend(['-t', input_cpt])
        
    run_gmx_command(grompp_args)
    
    mdrun_args = ['mdrun', '-v', '-deffnm', step_name]
    run_gmx_command(mdrun_args, cwd=output_dir)
    
    print(f"--- Equilibration step {step_num} finished ---")
    
    return {
        "gro": os.path.join(output_dir, f"{step_name}.gro"),
        "cpt": os.path.join(output_dir, f"{step_name}.cpt"),
        "edr": os.path.join(output_dir, f"{step_name}.edr"),
        "trr": os.path.join(output_dir, f"{step_name}.trr"),
    }
def run_equilibration_pipeline(topology_dir: str, output_dir: str):
    """
    Runs the full equilibration pipeline for a given topology.
    """
    print(f"--- Starting equilibration pipeline for {topology_dir} ---")
    
    os.makedirs(output_dir, exist_ok=True)
    
    base_name_list = [f.replace(".top", "") for f in os.listdir(topology_dir) if f.endswith(".top")]
    if not base_name_list:
        raise FileNotFoundError(f"No .top file found in {topology_dir}")
    base_name = base_name_list[0]
    
    # Initial files from minimization
    current_gro = os.path.join(topology_dir, "em.gro")
    topology_file = os.path.join(topology_dir, f"{base_name}.top")
    current_cpt = None
    
    # Use the energy-minimized structure as reference for restraints
    # (it has all atoms: protein + water + ions)
    reference_gro = current_gro  # Use em.gro as reference
    
    # MDP files for equilibration, in order
    mdp_files = [
        "config/2_equil.mdp",
        "config/3_equil.mdp",
        "config/4_equil.mdp",
        "config/5_equil.mdp",
        "config/6_equil.mdp",
        "config/7_prod_10ps_out.mdp",
    ]
    
    equilibration_results = []
    
    for i, mdp in enumerate(mdp_files):
        step_num = i + 2
        
        outputs = run_equilibration_step(
            step_num=step_num,
            mdp_file=os.path.abspath(mdp),
            input_gro=current_gro,
            input_top=topology_file,
            input_cpt=current_cpt,
            output_dir=output_dir,
            base_name=base_name,
            reference_gro=reference_gro  # Now points to em.gro
        )
        
        current_gro = outputs["gro"]
        current_cpt = outputs["cpt"]
        equilibration_results.append(outputs)
        
    print(f"--- Equilibration pipeline finished for {topology_dir} ---")
    return equilibration_results

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print("Usage: python -m scripts.protonation_pipeline.run_equilibration <topology_dir>")
        sys.exit(1)
        
    topology_dir = sys.argv[1]
    output_dir = f"{topology_dir}_equilibration"
    run_equilibration_pipeline(topology_dir, output_dir)