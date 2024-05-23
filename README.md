# PELE Constrained-Strain-Parallel (PELE_ConStraParl)

This repository is dedicated to preparing PELE adaptive simulations in MN5. It generates the file system needed to execute PELE ligand refinement simulations, starting from docked ligands.

In addition to regular PELE refinement simulations starting from docked ligands, this repository incorporates the following options:

- Parallel Execution: It can execute multiple ligand refinements in parallel against a single target.
- Predicted Binding Free Energy (BFE): It can compute the predicted PELE binding free energy for each ligand.

  $$BFE = \sum^N_i p_i E^b_i$$
  
  where N is the total number of accepted PELE steps, E^b_i the binding energy between the target and the ligand in the $i$ accepted PELE step, and p_i the probabilities obtained from a Boltzmann distribution using the N sampled complex (target-ligand) total energies:
  
  $$p_i=\frac{e^{E^t_i / KT}}{Q} \space \quad Q = \sum^{N}_{i}e^{E^t_i /KT}$$
  
- Constrained Refinement: Given a chemical constraint (currently, only hydrogen bonds are implemented), it can incorporate this constraint into the PELE refinement simulation of each ligand.
- Strain Energy Correction: It can adjust the binding free energies of the accepted PELE steps and the resulting PELE BFE by accounting for the ligand strain energy. To do so, short PELE simulations of the ligand with implicit solvent are performed to determine the ligand's minimal energy (using a temperature of 1500K and a single PELE epoch with 2,000 steps). For each accepted PELE step, the strain energy of the ligand is calculated by subtracting the minimal energy from the internal energy of the ligand. Finally, the corrected binding free energy for each accepted PELE step is obtained by adding the strain energy to the binding energy.

## Execution workflow:

* Input data:
  - **Docked ligands in a LIGS directory (ligands in individual PDB files without the target).**
  - **Prepared target protein (in PDB).**
  - (optional) Output table of the docking to filter ligands by properties in .csv format.
 
1) `prep_glides_to_PELE.py`
   **This script generates a COMPLEXES directory containing PDBs with the target and the ligands prepared for PELE.**
   (optional) It can create an HBList directory containing each complex's inter- and intra-molecular hydrogen bond interactions.
   (optional) It can produce an SDF file containing ligands that meet the specified property filters (requires the output table).

2) `generate_batch.py`
   **This script creates the execution files needed for each ligand to complete the PELE simulation and puts them into a new runs directory.**
   - runs_0_ligandX: SLURM script generating the subdirectory for the PELE simulations of ligandX.
   - batch_0: Batch file to execute all run_0_ligandX scripts in parallel.
   - runs_1_ligandX: SLURM script executing the PELE simulation of ligandX.
   - batch_1: Batch file to execute all run_1_ligandX scripts in parallel.
   - runs: Directory where all run and batch files will be stored.
   - results: Directory where all simulation output subdirectories and PELE configuration files will be stored.
  
3) `batch_0.sh`
   **This script generates the PELE configuration file (yaml PELE conf files) and the PELE simulation subdirectory within the results directory for each ligand.**
   It initiates PELE simulations subdirectories and halts PELE simulation.

4) `modify_conf_runs.py`
   **This script modifies the PELE configuration file (yaml PELE conf files) to include possible constraints (optional) and retrieve ligand internal energy for strain energy correction (optional) for each ligand.**
   Adds atom-atom constraints to pele.conf based on specified HB constraints.
   If strain correction is enabled, it incorporates the internal energy of the ligand into the report.
  
