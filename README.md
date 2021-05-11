# NFold
Process Machine Learning Histogram data to setup MELD folding simulations

## Instruction for running Meld simulation
04/30/2021 
PerezLab@UF

----------------------------------------
Required software:

OpenMM: https://github.com/openmm/openmm
Meld: https://github.com/maccallumlab/meld (Please make sure installing all required tools listed in README.md of Meld repository.)

-----------------------------------------

### step 1: Generate informational distograms, phi, psi from corresponding .npy prediction file.

```python analyze_distograms.py``` 

Description: Select high informational range of distograms, psi, phi from enormous amount of probability histogram predictions.
Input:  distogram.npy,  phi.npy,  psi.npy, sequence.fa                                      
Output: contacts.dat,  tight_contacts.dat,  phi.dat,  tight_phi.dat,  psi.dat,  tight_psi.dat

### step 2: Generate starting structure from Amber minimization.

```sbatch setup_from_random_cpu.sh``` 

Description: Before starting Meld simulation, we generate a starting system from given sequence and then minimize it with sander in Ambertools.
Input:  sequence.dat (contains only the second line of sequence.fa file)
Output: TEMPLATES/minimized.pdb

### step 3: Start Meld simulation with OpenMM.

```sbatch job_30_new.sh``` (contains 'python setup_aMeld.py' inside)

Description: Setup Meld simulation using ff14SBonlysc force field + restraints derived from output files in step 1.
Input:  setup_aMeld.py, sequence.dat, phi.dat, tight_phi.dat, psi.dat, tight_psi.dat, contact.dat, tight_contact.dat 
Output: Data/ remd.log (Data/ stores all simulation output, remd.log keep track of the simulation.)
