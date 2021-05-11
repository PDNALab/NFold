#!/bin/bash
#SBATCH --job-name=random            # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks=1      
#SBATCH --cpus-per-task=1            # Number of cores per MPI rank 
#SBATCH --nodes=1                    # Number of nodes
#SBATCH --ntasks-per-node=1          # How many tasks on each node
#SBATCH --ntasks-per-socket=1        # How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=800mb            
#SBATCH --distribution=cyclic:cyclic
#SBATCH --mem-per-cpu=800mb          # Memory per processor
#SBATCH --time=10:00:00              # Time limit hrs:min:sec
#SBATCH --output=random.log          # Standard output and error log


source ~/.load_OpenMM_cuda10         #load Ambertools, OpenMM+Meld
mkdir TEMPLATES                      #create directory to save minimized initial pdb file

cat<<EOF>setup_random.py
#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld import system
from meld import parse

def setup_system():
    # load the sequence
    sequence = parse.get_sequence_from_AA1(filename='sequence.dat')
    n_res = len(sequence.split())
    
    # build the system
    p = system.ProteinMoleculeFromSequence(sequence)
    b = system.SystemBuilder(forcefield="ff14sbside")
    s = b.build_system_from_molecules([p])
setup_system()
EOF

python  setup_random.py > tleap.in
tleap -f tleap.in

cat<<EOF>minimize.in
Stage 1 - minimisation of 1sr protein
 &cntrl
  imin=1, maxcyc=1000, ncyc=500,
  cut=999., rgbmax=999.,igb=8, ntb=0,
  ntpr=100,
 /
EOF
source ~/.load_Amber
sander -O \
    -i minimize.in \
    -o minimize.out \
    -c system.mdcrd \
    -ref system.mdcdr \
    -r eq0.rst \
    -p system.top \
    -e eq0.ene \
    -x eq0.netcdf

cpptraj system.top<<EOF>&err
trajin eq0.rst
trajout TEMPLATES/minimized.pdb pdb
go
EOF

