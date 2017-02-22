#!/bin/bash
#SBATCH --job-name=as_2Y95
#SBATCH --time=72:00:00
#SBATCH --ntasks=18
#SBATCH --ntasks-per-node=18
#SBATCH -p brooks

# Load Modules 
source ~/.bashrc
export PATH=/home/afrank/local_software/bin/:$PATH
cd /home/afrankz/testbed/testWilliam/calcs

rna=2Y95
./scripts/assigner.R --output=assignments/assigned_shifts_ramsey_${rna} --iterations=10000 --freq_output=1000 --scale=0.75 --parallel --nprocessors=18 data/ramsey_${rna}.txt data/observed_shifts_${rna}.txt data/ramsey_accuracy_resname_nucleus.txt


