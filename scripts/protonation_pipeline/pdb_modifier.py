import os
import subprocess
from typing import Dict, Tuple

def run_propka(pdb_file: str) -> str:
    """Runs PROPKA on a PDB file and returns the path to the .pka file."""
    print(f"Running PROPKA on {pdb_file}...")
    
    # Check if propka is installed
    try:
        subprocess.run(["propka3", "--version"], check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        raise RuntimeError("PROPKA is not installed or not in PATH. Please install it to proceed.")

    propka_command = ["propka3", pdb_file]
    subprocess.run(propka_command, check=True, capture_output=True)
    
    pka_file = os.path.basename(pdb_file).replace(".pdb", ".pka")
    if not os.path.exists(pka_file):
        # propka 3.1 puts it in the same dir, but let's be safe
        pka_file = pdb_file.replace(".pdb", ".pka")

    if not os.path.exists(pka_file):
        raise FileNotFoundError(f"PROPKA output file not found: {pka_file}")

    print(f"PROPKA finished. Output written to {pka_file}")
    return pka_file

def parse_pka(pka_file: str) -> Dict[Tuple[str, int], float]:
    """Parses a .pka file and returns a dictionary of pKa values."""
    pka_data = {}
    with open(pka_file, "r") as f:
        for line in f:
            if not line.startswith(" "):
                continue
            parts = line.split()
            if len(parts) >= 4:
                try:
                    res_name = parts[0]
                    res_id = int(parts[1])
                    chain_id = parts[2]
                    pka = float(parts[3])
                    pka_data[(res_name, res_id, chain_id)] = pka
                except (ValueError, IndexError):
                    continue
    return pka_data

def modify_pdb_for_ph(pdb_file: str, pka_data: Dict[Tuple[str, int], float], target_ph: float) -> str:
    """Creates a new PDB file with residue names modified for a target pH."""
    
    output_pdb = os.path.basename(pdb_file).replace(".pdb", f"_pH{int(target_ph)}.pdb")
    
    # Create a mapping of (res_id, chain_id) to new residue name
    residue_map = {}
    for (res_name, res_id, chain_id), pka in pka_data.items():
        new_res_name = res_name
        if res_name in ["ASP", "GLU"] and target_ph < pka:
            new_res_name = "ASH" if res_name == "ASP" else "GLH"
        elif res_name == "HIS":
            if target_ph < pka:
                new_res_name = "HIP" # Doubly protonated
            else:
                # Neutral HIS can be HID or HIE. Defaulting to HID.
                # A more sophisticated approach would check H-bonding.
                new_res_name = "HID"
        
        if new_res_name != res_name:
            residue_map[(str(res_id), chain_id)] = new_res_name

    with open(pdb_file, "r") as infile, open(output_pdb, "w") as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_id = line[22:26].strip()
                chain_id = line[21]
                if (res_id, chain_id) in residue_map:
                    original_res_name = line[17:20]
                    new_res_name = residue_map[(res_id, chain_id)]
                    line = line[:17] + new_res_name + line[20:]
            outfile.write(line)
            
    print(f"Created modified PDB file for pH {target_ph}: {output_pdb}")
    return output_pdb
