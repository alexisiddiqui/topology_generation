import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np

def extract_energy_terms(edr_file: str, terms: list, output_dir: str):
    """Extracts energy terms from an EDR file using gmx energy."""
    print(f"Extracting {', '.join(terms)} from {edr_file}...")
    
    output_xvg = os.path.join(output_dir, f"{os.path.basename(edr_file)}.xvg")
    
    # Prepare input for gmx energy (term numbers followed by 0)
    gmx_input = b"\n".join([term.encode() for term in terms]) + b"\n0\n"
    
    gmx_energy = subprocess.run(
        ['gmx', 'energy', '-f', edr_file, '-o', output_xvg],
        input=gmx_input,
        capture_output=True,
        text=True
    )
    
    if gmx_energy.returncode != 0:
        print(f"Error running gmx energy: {gmx_energy.stderr}")
        return None
        
    return output_xvg

def calculate_rmsd(tpr_file: str, trr_file: str, output_dir: str):
    """Calculates RMSD using gmx rms."""
    print(f"Calculating RMSD for {trr_file}...")
    
    output_xvg = os.path.join(output_dir, f"rmsd_{os.path.basename(trr_file)}.xvg")
    
    # Input for gmx rms: Backbone for fitting, Backbone for RMSD calc
    gmx_input = b"Backbone\nBackbone\n"
    
    gmx_rms = subprocess.run(
        ['gmx', 'rms', '-s', tpr_file, '-f', trr_file, '-o', output_xvg, '-tu', 'ps'],
        input=gmx_input,
        capture_output=True,
        text=True
    )
    
    if gmx_rms.returncode != 0:
        print(f"Error running gmx rms: {gmx_rms.stderr}")
        return None
        
    return output_xvg

def plot_xvg(xvg_files: list, title: str, xlabel: str, ylabel: str, output_file: str):
    """Plots data from one or more XVG files."""
    plt.figure(figsize=(10, 6))
    
    for xvg_file in xvg_files:
        data = np.loadtxt(xvg_file, comments=["@", "#"])
        label = os.path.basename(xvg_file).split('_')[0]
        plt.plot(data[:, 0], data[:, 1], label=label)
        
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.grid(True)
    plt.savefig(output_file, dpi=300)
    plt.close()

def analyze_equilibration(equilibration_dir: str, output_dir: str):
    """
    Analyzes the equilibration results.
    """
    print(f"--- Analyzing equilibration data in {equilibration_dir} ---")
    os.makedirs(output_dir, exist_ok=True)
    
    edr_files = sorted([os.path.join(equilibration_dir, f) for f in os.listdir(equilibration_dir) if f.endswith(".edr")])
    trr_files = sorted([os.path.join(equilibration_dir, f) for f in os.listdir(equilibration_dir) if f.endswith(".trr")])
    tpr_files = sorted([os.path.join(equilibration_dir, f) for f in os.listdir(equilibration_dir) if f.endswith(".tpr")])

    # --- Plot Temperature ---
    temp_xvgs = [extract_energy_terms(f, ["Temperature"], output_dir) for f in edr_files]
    plot_xvg(temp_xvgs, "Temperature vs. Time", "Time (ps)", "Temperature (K)", os.path.join(output_dir, "temperature.png"))

    # --- Plot Pressure ---
    pressure_xvgs = [extract_energy_terms(f, ["Pressure"], output_dir) for f in edr_files]
    plot_xvg(pressure_xvgs, "Pressure vs. Time", "Time (ps)", "Pressure (bar)", os.path.join(output_dir, "pressure.png"))

    # --- Plot Density ---
    density_xvgs = [extract_energy_terms(f, ["Density"], output_dir) for f in edr_files]
    plot_xvg(density_xvgs, "Density vs. Time", "Time (ps)", "Density (kg/m^3)", os.path.join(output_dir, "density.png"))

    # --- Plot RMSD ---
    rmsd_xvgs = [calculate_rmsd(tpr, trr, output_dir) for tpr, trr in zip(tpr_files, trr_files)]
    plot_xvg(rmsd_xvgs, "RMSD vs. Time", "Time (ps)", "RMSD (nm)", os.path.join(output_dir, "rmsd.png"))
    
    print(f"--- Analysis complete. Plots are in {output_dir} ---")

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print("Usage: python -m scripts.protonation_pipeline.analyze_equilibration <equilibration_dir>")
        sys.exit(1)
        
    equilibration_dir = sys.argv[1]
    output_dir = f"{equilibration_dir}_analysis"
    analyze_equilibration(equilibration_dir, output_dir)
