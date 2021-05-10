#!/bin/bash
#SBATCH --job-name=1dj7B00               # Job name
#SBATCH --mail-type=END,FAIL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=30      
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-gpu=1
#SBATCH --partition=gpu
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=800mb              # Memory per processor
#SBATCH --time=60:00:00                  # Time limit hrs:min:sec
#SBATCH --output=meld.log                # Standard output and error log

 


source ~/.load_OpenMM_cuda10             #load OpenMM+Meld

[[ -d Data ]] || python setup_aMeld.py   #check if there is already a Data/, we are continuing a killed simulation, otherwise start new setup_aMeld.py simulation.

if [ -e remd.log ]; then                 #First check if there is a remd.log file, we are continuing a killed simulation
    prepare_restart --prepare-run        #so we need to prepare_restart.
      fi

srun --mpi=pmix  launch_remd --debug     #restart remd simulation
