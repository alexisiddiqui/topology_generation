# Protein Topology Generation and Equilibration Pipeline

This project provides a pipeline for generating GROMACS topologies for proteins at different pH values, followed by an equilibration simulation and analysis.

## Overview

The pipeline automates the following steps:
1.  **PDB Download**: Downloads a protein structure from the RCSB PDB database.
2.  **Protonation State Prediction**: Uses PROPKA to predict the pKa values of titratable residues.
3.  **Topology Generation**: Creates GROMACS topology files for the protein at specified pH values (e.g., pH 7 and pH 4). This includes:
    *   Modifying residue names based on their predicted protonation state.
    *   Generating the topology using `pdb2gmx`.
    *   Creating a simulation box, solvating it, and adding ions.
    *   Performing an energy minimization.
4.  **Equilibration**: Runs a multi-step equilibration simulation to relax the system.
5.  **Analysis**: Analyzes the equilibration trajectory to check for stability, calculating and plotting:
    *   Temperature
    *   Pressure
    *   Density
    *   RMSD

## Requirements

### Python Packages

You can install the required Python packages using `pip`. It's recommended to do this in a virtual environment.

```bash
pip install mdanalysis matplotlib numpy requests propka
```

### Command-Line Tools

*   **GROMACS**: The pipeline scripts call GROMACS executables (`gmx` or `gmx_mpi`) directly using `subprocess`. GROMACS must be installed and sourced in your environment. See the installation guide below.
*   **PROPKA**: This tool is used for predicting pKa values. The `pip install` command above will install the `propka` command line tool.

## GROMACS Installation Guide (with CUDA via Conda)

This guide details how to install a GPU-accelerated build of GROMACS in a dedicated Conda environment. This is the recommended method for this project.

**System Information:**
*   **GPU:** NVIDIA GeForce RTX 3090
*   **CUDA Driver:** 12.2

**Installation Steps:**

1.  **Create and Activate a New Conda Environment**

    This creates an isolated environment named `gromacs-cuda` to avoid conflicts with system packages.

    ```bash
    conda create -n gromacs-cuda -c conda-forge
    conda activate gromacs-cuda
    ```

2.  **Install GROMACS with CUDA**

    This command installs GROMACS and a compatible CUDA toolkit from the `conda-forge` channel. Conda will handle the specific CUDA versioning for you.

    ```bash
    conda install -c conda-forge gromacs cudatoolkit
    ```

3.  **Verify the GROMACS Installation**

    Check that the newly installed GROMACS is being used and that it was built with CUDA support.

    ```bash
    which gmx 
    # Should point to the anaconda directory
    
    gmx --version
    # Look for "CUDA support:      enabled" in the output
    ```

## GROMACS Installation Guide (with CUDA via venv/source)

This guide details how to install a GPU-accelerated build of GROMACS from source, which allows for more control and is suitable for use with a standard Python `venv`.

**Prerequisites:**
*   A C++ compiler (e.g., `gcc`/`g++`)
*   CMake
*   NVIDIA CUDA Toolkit (installed on your system)

**Installation Steps:**

1.  **Download and Extract GROMACS Source Code**

    Download the desired version of GROMACS from the official website (e.g., version 2023.3).

    ```bash
    wget https://ftp.gromacs.org/gromacs/gromacs-2023.3.tar.gz
    tar -xzvf gromacs-2023.3.tar.gz
    cd gromacs-2023.3
    ```

2.  **Create a Build Directory**

    It's best practice to build GROMACS in a separate directory.

    ```bash
    mkdir build
    cd build
    ```

3.  **Configure the Build with CMake**

    This command configures the build to use CUDA and sets the installation directory. Replace `/path/to/your/gromacs_install` with your desired installation location.

    ```bash
    cmake .. -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=/path/to/your/gromacs_install
    ```

4.  **Compile and Install GROMACS**

    This will compile and install GROMACS to the prefix you specified.

    ```bash
    make
    make check # Optional, but recommended
    make install
    ```

5.  **Source GROMACS in Your Environment**

    To use the compiled GROMACS, you need to source the `GMXRC` script from your installation directory. You should add this line to your shell's startup file (e.g., `~/.bashrc`) to make it permanent.

    ```bash
    source /path/to/your/gromacs_install/bin/GMXRC
    ```

6.  **Verify the Installation**

    Activate your Python `venv`, source the `GMXRC` file if you haven't already, and verify the installation.

    ```bash
    which gmx
    # Should point to your custom installation path
    
    gmx --version
    # Look for "CUDA support:      enabled" in the output
    ```

## How to Run the Pipeline

The entire pipeline can be run with a single command:

```bash
python run_pipeline.py <pdb_id>
```

Replace `<pdb_id>` with the 4-character PDB ID of the protein you want to process (e.g., `2L39`).

The pipeline will create several directories:
*   `input_structures/<pdb_id>/`: Contains the downloaded PDB file.
*   `top/<pdb_id>_pH<ph>/`: Contains the generated topology files.
*   `top/<pdb_id>_pH<ph>_equilibration/`: Contains the files from the equilibration simulation.
*   `top/<pdb_id>_pH<ph>_equilibration_analysis/`: Contains the analysis plots.

## Pipeline Scripts

The `scripts/protonation_pipeline/` directory contains the individual modules that make up the pipeline:

*   `download_pdb.py`: Downloads a PDB file from the RCSB database.
*   `pdb_modifier.py`: Runs PROPKA to determine protonation states of residues at a given pH and modifies the PDB file accordingly.
*   `generate_topology.py`: Uses the modified PDB file to generate a GROMACS topology, create a simulation box, add solvent and ions, and run an energy minimization. It calls GROMACS command-line tools directly.
*   `run_equilibration.py`: Executes a series of equilibration steps (NVT and NPT) to prepare the system for a production simulation. It calls GROMACS command-line tools directly.
*   `analyze_equilibration.py`: Analyzes the results of the equilibration, plotting key metrics like temperature, pressure, density, and RMSD to ensure the system is stable.