# PELE Constrained-Strain-Parallel (PELE_ConStraParl)

This repository is dedicated to preparing PELE adaptive simulations in MN5. It generates the file system needed to execute PELE ligand refinement simulations, starting from docked ligands.

## Repository features:

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
 
1) `prep_glides_to_PELE.py`<br>
   **This script generates a COMPLEXES directory containing PDBs with the target and the ligands prepared for PELE.** <br>
   (optional) It can create an HBlist directory containing target and target-ligand hydrogen bond interactions for each ligand. <br>
   (optional) It can produce an SDF file containing ligands that meet the specified property filters (requires the output table).

3) `generate_batch.py`<br>
   **This script creates the execution files needed for each ligand to complete the PELE simulation and puts them into a new runs directory.** <br>
   - runs_0_ligandX: SLURM script generating the subdirectory for the PELE simulations of ligandX.
   - batch_0: Batch file to execute all run_0_ligandX scripts in parallel.
   - runs_1_ligandX: SLURM script executing the PELE rescoring simulation of ligandX.
   - batch_1: Batch file to execute all run_1_ligandX scripts in parallel.
   - runs_2_ligandX (only if strain): SLURM script executing the PELE simulation of the ligand alone in implicit solvent.
   - batch_2 (only if strain): Batch file to execute all run_2_ligandX scripts in parallel.
   - runs_3_ligandX (only if strain): SLURM script applying strain corrections.
   - batch_3 (only if strain): Batch file to execute all run_3_ligandX scripts in parallel.
   - runs: Directory where all run and batch files will be stored.
   - results: Directory where all simulation output subdirectories and PELE configuration files will be stored.
  
4) `batch_0.sh`<br>
   **This script generates the PELE configuration file (yaml PELE conf files) and the PELE simulation subdirectory within the results directory for each ligand.** <br>
   It initiates PELE simulations subdirectories and halts PELE simulation.

5) `modify_conf_runs.py` <br>
   **This script modifies the PELE configuration file (yaml PELE conf files) to include possible constraints (optional) and retrieve ligand internal energy for strain energy correction (optional) for each ligand.** <br>
   Adds atom-atom constraints to pele.conf based on specified HB constraints. You can add more than one constraint and only the ones found in each ligand would by applied.
   If strain correction is enabled, it incorporates the internal energy of the ligand into the report.

6) `batch_1.sh` <br>
   **This script will start all PELE rescoring simulations** <br>

7) `batch_2.sh` <br>
   **If strain correction is enabled, this script will run short PELE simulations, one for each compound (to get its minimal energy)** <br>
   <ins>This script can be run simultaneously with batch_1.sh script!</ins>

8) `batch_3.sh` <br>
   **If strain correction is enabled, this script will apply it to all binding free energies of the rescoring simulation and the resulting PELE BFE** <br>
9) `pele_analysis.py` <br>
  **This script generates a CSV file containing PELE BFE-related parameters for each ligand. It also retrieves a representative pose of the simlated ligand (the accepted PELE step with the lowest binding energy) with the original connectivity.** <br>
(Optional) If a CSV file with docking parameters is provided, it will include additional metrics in the final CSV file.
