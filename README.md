# particle_in_mag
Reproduction of geant4 magnetic fields.
Tue file `particle_in_mag` contains the source code to reproduce the siumulation of a particle traveling through air with magnetic field. It receives the initial condition of the particle and outputs the trajectory [(x, y, z, px, py, pz)] of this particle after N steps.

## Geant4 simulations

`geant4.py` contains the core simulation code for modeling particle interactions within magnetic fields in air using the Geant4 toolkit. The Geant4 implementation is presented in `muons_and_matter_cpp`, based on the library MuonsAndMatter for simulating muons travelling through space using Geant4.

## Trajectories comparison

`comparte_trajectories.py` contains a script for running both simulations and comparing the trajectories. 3 different magnetic fields are presented, Toy magnetic field (arbitrary implementation for testing), Uniform, where a fixed B is set for the entire space and Map, where the magnetic field is taken from some pickle file.
