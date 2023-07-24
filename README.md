CC LLPS Framework
==============================
Custom analysis code and MD parameter files related to the manuscript titled "Coiled-coil domains are sufficient to drive liquid-liquid phase separation of proteins in molecular simulation". Authors: Dominique A. Ramirez, Loren E. Hough, Michael R. Shirts
The code in this repository are under the MIT License.


## Description of provided files

### Custom\_scripts
This folder contains the scripts to do the multimerization counting analysis, which was used in Figure S4. In order to use `multimerization_analysis.py`, a GROMACS MD trajectory file (.xtc) and an associated topology file (`.pdb`) are needed. You also need a "Protein Building" design file, which specifies how the CC proteins in the simulation are organized, and this file is described in more detail in the `./ProteinBuilding/proteinBuilding_README.md` file within this repository. This code calculates the numbers and types of multimers in each specified frame of the trajectory, and then produces an averaged distribution of multimers and a timeseries of each type of multimer. These scripts were written specifically for analysis of simulation trajectories for this paper.

Other analyses like density profile and molecular cluster analyses are tools provided by GROMACS. 


### MD\_parameter\_files
This folder contains the `.mdp` files for the single molecule and slab simulation procedures.


### ProteinBuilding
This folder contains the scripts to do the actual protein building according to the `PeptideBuilder` strategy described in Supporting Material, Section S.VIII. Instructions and an example are provided to demonstrate how the code works.

As described in the paper, once a protein is made, you can then run single molecule simulations on it using the appropriate `.mdp` files (described above) and then pack configurations into a box of desired density and size using `Packmol v. 20.3.2`.


### Topologies
This folder contains all of the topology files (`.itp` and `.top`) for all of the proteins described in this manuscript. Each protein pair is given its own folder. `.top` files are predesigned for slab simulations.


### Snapshots
This folder contains snapshots of slab images that appear in the manuscript. Image files (.png) and PyMOL sessions (.pse) are provided.


## Guidelines for reproducing data
Due to the large size of individual simulation trajectories (&#076; 3 Gb), trajectories are not provided. You can follow these guidelines to reproduce trajectory data.
1. Write a protein `design_file.csv` (see README in the `ProteinBuilding` directory) for whichever protein you want.
2. Run a single molecule simulation at the desired slab temperature using the provided `.mdp` files. 
    - takes less than 1 hour using 1 thread on a 2019 Macbook Pro with a 2 Quad-core Intel CPU
3. Randomly select 5 configurations from the equilibrated portion of the single molecule simulation (using whatever random selection method you'd like).
4. Pack enough copies of all 5 configurations to reach the desired coil segment density in the box.
    - make the total number of *coil segments* in the box equal to 450 to be consistent with experiments in this manuscript.
5. Run the slab protocol (EM--NPT--NVT--Production MD steps) using the provided `.mdp` files.
    - EM--NPT--NVT steps take less than 8 hours on a high performance computer using 4 processes
    - the slab production step takes, on average, less than 7 days on a high performance computer using 24 processes depending on available resources
6. Use `gmx_density` to determine the equilibrated portion of the trajectory, then to do profile analysis, and `gmx_clustsize` to calculate cluster size information
7. Use `gmx_msd` to calculate the mean squared displacement on all proteins in the system. Bootstrapping can then be done using the MSD data to estimate an effective diffusion coefficient from the early linear portion of the MSD vs. lag-time data.
    - `gmx_msd` requires extensive memory usage, but finishes quickly. For a 20 microsecond trajectory, this analysis takes less than 10 minutes on a high performance computer with at least 30 GB of allocated RAM


## Software requirements and details

### Molecular simulation software
We used GROMACS v2022.1 to run all molecular dynamics simulations, both single molecule and slab, for this study. The specific details of compilation and types of computers used are described in the manuscript both in Methods and Supplemental Methods. 

You may download GROMACS v2022.1 [at this link](https://manual.gromacs.org/2022.1/download.html), and you can reference the documentation for this specific version [at this link](https://manual.gromacs.org/2022.1/index.html). GROMACS compilation is hardware specific and so we cannot provide precompiled binaries of GROMACS. Installation is quick (within 2 hours) and is easy to do following the instructions provided in the documentation.

We used following compilers were used to compile GROMACS:
- MacOS (Clang 12.0.0)
- Ubuntu (GCC 11.3.0)
- RHEL (GCC 11.2.0; OpenMPI 4.1.1)


### Requirements to run custom scripts
Custom scripts -- which include the custom analysis scripts (see `Custom_scripts`, above) and the scripts used to make proteins (see `ProteinBuilding`, above) -- use python v3.9. Other versions of python have not been tested. The following dependencies are also needed:
- matplotlib (3.7.1)
- numpy (1.21.5)
- mdtraj (1.9.7)
- ProDy (2.3.1)
- PeptideBuilder (1.1.0)
- Biopython (1.80)
- pandas (1.4.4)
- scipy (1.9.1)

These scripts were run on the following operating systems:
- MacOS 12.6.3
- Ubuntu 22.04.2 LTS
- RHEL 7.9 (codename Maipo)


### Installation guide
Running the custom scripts requires installing the above listed dependencies (and their own subsequent dependencies). Installing GROMACS requires following the instructions provided in its documentation (see link above). GROMACS installation is the limiting step, and should take less than 2 hours on a normal computer. If you have access to a high performance computer, GROMACS might already be installed.



### Copyright

Copyright (c) 2023, Dominique Ramirez
