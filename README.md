# MPhil_Project_Chromatin: Phase transitions of our genome with multiscale coarse-grained simulations, supplementary code

This project investigates the phase separation of chromatin. It employes the Collepardo group's minimal chromatin model (https://github.com/CollepardoLab/CollepardoLab_Chromatin_Model.git), 
running coexistence simulations with LLAMPS (https://github.com/lammps/lammps.git) and calculating the critical point and phase diagram using my own code.

## Simulation setup files
The files for the simulation setup are at the LAMMPS_simulation_setup folder. It contains:
- create_slab.in:                   Initializes and condenses nucleosome arrays.
- run_slab_mixture.in:              Equilibrates the system.
- data_files.in:                    Reads pre-relaxed nucleosome arrays obtained from replica exchange simulations
                                    (it was shared by the Collepardo lab with due to high computational costs and the scope of the project)
- DNA_sequence.txt:                 Contains linker DNA and histone-DNA sequences.
- NAFlex_params.txt:                Stores helical equilibrium parameters and stiffness matrix of the chromatin model.
- chromatin.so:                     Plugin implementing the Collepardo chromatin model (compiled during LAMMPS build).
- create_slub.sub, run_slab.sub:    Submission scripts for running LAMMPS jobs on the Icelake node of CSD3.
- create_simulation.sh:             Shell script to initialize simulations at user-defined salt concentrations.

## Creating the environment
The environment was step following the instructions below. The echod $LD_LIBRARY_PATH was included in the .sub submission scripts

git clone https://github.com/lammps/lammps.git
cd lammps
git checkout tags/stable_29Sep2021 -b stable
mkdir build && cd build
cmake -D BUILD_SHARED_LIBS=yes -D BUILD_TOOLS=yes -D PKG_ASPHERE=yes -D PKG_RIGID=yes -D PKG_MOLECULE=yes -D PKG_PLUGIN=yes -D PKG_EXTRA-PAIR=yes -D PKG_COLVARS=yes -D PKG_EXTRA-FIX=yes -DCMAKE_CXX_FLAGS='-O3 -march=native -fno-math-errno -ffast-math' ../cmake
make install
git clone https://github.com/CollepardoLab/CollepardoLab_Chromatin_Model.git
cd CollepardoLab_Chromatin_Model/
ls
git checkout plugin-dev
mkdir build && cd build
cmake -DLAMMPS_HEADER_DIR=/rds/user/jip29/rds-t2-cs166-GfNGl1zzWsk/NACHO/lammps/src/ ../../CollepardoLab_Chromatin_Model/
make install
echo $LD_LIBRARY_PATH

Other dependences include numpy, scipy, matplotlib, pandas

## Processing the simulation
The code calc_phase_diagram.py processes the simulation output trjaectories (.dump) and calculates the phase diagram of the chromatin system.
- Input: name_dir (The name of the directory containing the output .dump files, with the value of the simulation salt
                    concentraion appended in their name. For example, cores_0.074.dump)
- Outputs:
  - density_profile_and_fit: Folder containing the fitted density profiles in a .png plot and the density values in a .txt
  - densities.csv: Csv file containing the salt concentration of the simulation, the extracted low and extracted high densities of the simulaiton in rows
  - phase_diagram.png:  The phase diagram of the simulatin
  - 

