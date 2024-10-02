# DCVGCMD

This is an test case for the DCVGCMD method implement in Yuan et al...

It shoud contain three folders:

src, data, test-pim-1, and a readme file. 

## Prerequisites

To perform the simulation, you will need

Python3, and GROMACS (2023.3 or similar)

## Running

To run the DCVGCMD, go to the directory test-pim-1,

```
cd test-pim-1
```

and source the run.sh file

```
source run.sh
```

Please be careful with the settings and parameters. Make sure the gromacs prefix is the same as your local environment and check if it works with -ntmpi. 

## Useful information

### Output: 

The folder genereated with numbers are the simulations from the DCVGCMD setup. There are several files that you should find useful other than the standard gromacs files, e.g., trajectories (.xtc), energies (.edr). 

Two info files: info_inlet and info_outlet contains the atom addition/deletion of the inlet and outlet CVs, respectively. When calculating the permeability, info_outlet will be useful. 

The posres.gro file is the file used for position restraining in GROMACS. 

