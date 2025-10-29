#!/bin/bash
#SBATCH --account=brave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared
#SBATCH --time=24:00:00
#SBATCH --job-name=min
#SBATCH --out=slurm%j.out
#SBATCH --mem-per-cpu=10GB


module purge
module load mamba cuda/12.2
source /projects/robustmicrob/jlaw/tools/amber24/amber.sh

pmemd -O -i ../../script_25/min1.in -o min1.out -p complex_solvated.prmtop -c complex_solvated.inpcrd -r min1.rst -inf min1.mdinfo -ref complex_solvated.inpcrd
pmemd -O -i ../../script_25/min2.in -o min2.out -p complex_solvated.prmtop -c min1.rst -r min2.rst -inf min2.mdinfo -ref min1.rst


cp complex_solvated.prmtop min2.rst system_hmass.top trial_1_boltz
cp complex_solvated.prmtop min2.rst system_hmass.top trial_2_boltz
cp complex_solvated.prmtop min2.rst system_hmass.top trial_3_boltz
