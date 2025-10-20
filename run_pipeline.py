import sys
from scripts.protonation_pipeline.download_pdb import download_pdb
from scripts.protonation_pipeline.generate_topology import generate_topology_for_ph
from scripts.protonation_pipeline.run_equilibration import run_equilibration_pipeline
from scripts.protonation_pipeline.analyze_equilibration import analyze_equilibration

def main():
    """
    Main function to run the topology generation pipeline.
    """
    if len(sys.argv) != 2:
        print("Usage: python run_pipeline.py <pdb_id>")
        sys.exit(1)
        
    pdb_id = sys.argv[1]
    
    print(f"Starting topology generation pipeline for PDB ID: {pdb_id}")
    
    # 1. Download PDB file
    pdb_file = download_pdb(pdb_id, output_dir=f"input_structures/{pdb_id}")
    
    if not pdb_file:
        print(f"Could not download PDB file for {pdb_id}. Exiting.")
        sys.exit(1)
        
    # 2. Generate topologies for pH 7 and pH 4
    for ph in [7.0, 4.0]:
        try:
            topology_dir = f"top/{pdb_id}_pH{int(ph)}"
            generate_topology_for_ph(pdb_file, ph=ph, output_dir=topology_dir)
            
            equilibration_dir = f"{topology_dir}_equilibration"
            run_equilibration_pipeline(topology_dir, equilibration_dir)
            
            analysis_dir = f"{equilibration_dir}_analysis"
            analyze_equilibration(equilibration_dir, analysis_dir)

        except Exception as e:
            print(f"An error occurred during pipeline for pH {ph}: {e}")
            continue # Continue to the next pH value
        
    print("\nPipeline finished successfully!")
    print(f"Results for pH 7 are in: top/{pdb_id}_pH7_equilibration_analysis")
    print(f"Results for pH 4 are in: top/{pdb_id}_pH4_equilibration_analysis")

if __name__ == "__main__":
    main()

