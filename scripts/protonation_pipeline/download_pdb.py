import os
import requests

def download_pdb(pdb_id, output_dir="input_structures"):
    """
    Downloads a PDB file from the RCSB PDB database.

    Args:
        pdb_id (str): The 4-character PDB ID.
        output_dir (str): The directory to save the PDB file.

    Returns:
        str: The path to the downloaded PDB file.
    """
    pdb_id = pdb_id.upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    os.makedirs(output_dir, exist_ok=True)
    
    output_path = os.path.join(output_dir, f"{pdb_id}.pdb")

    if os.path.exists(output_path):
        print(f"PDB file for {pdb_id} already exists at {output_path}")
        return output_path

    print(f"Downloading PDB file for {pdb_id} from {url}...")
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        with open(output_path, "w") as f:
            f.write(response.text)
            
        print(f"Successfully downloaded {pdb_id}.pdb to {output_path}")
        return output_path

    except requests.exceptions.RequestException as e:
        print(f"Error downloading PDB file: {e}")
        return None

if __name__ == '__main__':
    download_pdb("2L39")
