# PELE Constrained-Strain-Parallel (PELE_ConStraParl)

This repository is dedicated to preparing PELE adaptive simulations in MN5. It generates the file system needed to execute PELE ligand refinement simulations, starting from docked ligands.

In addition to regular PELE refinement simulations starting from docked ligands, this repository incorporates the following options:

- Parallel Execution: It can execute multiple ligand refinements in parallel against a single target.
- Predicted Binding Free Energy (BFE): It can compute the predicted PELE binding free energy for each ligand.

  $$BFE = \sum^N_i p_i$$
  
  where N is the total number of accepted PELE steps, E^b_i the binding energy between the target and the ligand in the $i$ accepted PELE step, and p_i the probabilities obtained from a Boltzmann distribution using the N sampled complex (target-ligand) total energies:
  
  $$p_i=\frac{e^{E^t_i / KT}}{Q} \space \quad Q = \sum^{N}_{i}e^{E^t_i /KT}$$
  
- Constrained Refinement: Given a chemical constraint (currently, only hydrogen bonds are implemented), it can incorporate this constraint into the PELE refinement simulation of each ligand.
- Strain Energy Correction: It can adjust the binding free energies of the accepted PELE steps and the resulting PELE BFE by accounting for the ligand strain energy. To do so, short PELE simulations of the ligand with implicit solvent are performed to determine the ligand's minimal energy (using a temperature of 1500K and a single PELE epoch with 2,000 steps). For each accepted PELE step, the strain energy of the ligand is calculated by subtracting the minimal energy from the internal energy of the ligand. Finally, the corrected binding free energy for each accepted PELE step is obtained by adding the strain energy to the binding energy.

## Execution workflow:

