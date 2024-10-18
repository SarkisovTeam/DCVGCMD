# Dual Control Volume Grand Canonical Molecular Dyanmics

This repository contains a test case for the Dual Control Volume Grand Canonical Molecular Dynamics (DCVGCMD) method, as implemented by Yuan et al. The method allows for molecular dynamics simulations with controlled exchange of atoms between control volumes, simulating processes like gas permeating through membranes.

## Repository Structure
The repository is organized into the following key folders:

- **`src/`**: Contains Python scripts and tools for running the simulation.
- **`data/`**: Includes forcefield files and molecular data necessary for simulations (e.g., PIM-1, CO2, and N2).
- **`test-pim-1/`**: Example test case setup for the simulation, including all necessary files to run a sample simulation.

## Prerequisites

Before running the DCVGCMD simulation, make sure you have the following software installed:

- **Python 3.x**
- **GROMACS** (version 2023.3 or similar)

Ensure that your GROMACS environment is correctly configured and accessible from your terminal.


## How to Run the Simulation

To run the DCVGCMD simulation using the provided test case:

1. Navigate to the `test-pim-1` directory:

   ```bash
   cd test-pim-1
   ```

2. Source the `run.sh` file to set up the environment and initiate the simulation:
  
   ```bash
   source run.sh
   ```

### Important Notes:
-  **GROMACS Prefix**: Ensure that the GROMACS prefix in the script matches your local GROMACS installation.
-  **Number of Threads**: Check if the -ntmpi option works with your GROMACS version to manage the number of threads efficiently.


## Useful information

### Source Code: `src/`

-  **replace-cv-5s-v5.py** This script handles the exchange of control volumes during the simulation.
You can view usage options by running:



```bash
python replace-cv-5s-v5.py -h 
```

An example of how to use the script is provided in the `test-pim-1` folder.



### Output: 

During the simulation, the following outputs will be generated:

- **Simulation Folders**: Folders with numerical names corresponding to different simulation runs (starting from 0).
- **Trajectory Files**: .xtc files containing simulation trajectories.
- **Energy Files**: .edr files storing energy data from the simulation.
- **Info Files**:
  - `info_inlet`: Contains information about atom additions and deletions in the inlet control volume.
  - `info_outlet`: Similar to `info_inlet`, but for the outlet control volume. This file is especially useful when calculating permeability.
- **Position Restraint File**: `posres.gro` – This file is used for position restraining during GROMACS simulations.



### Data Files

The `data/` folder contains the following key files:

- **Forcefield files**: For PIM-1, CO2, and N2 molecules.
- **Reservoir Example**: A bulk reservoir example with a 15:85 CO2 mixture at 1 bar, ready for use in the simulation.


## Contact

If you encounter any issues or have questions, please contact:



**Tianmu (Tim) Yuan**

Department of Chemical Engineering

The University of Manchester

tianmu.yuan(at)manchester.ac.uk
