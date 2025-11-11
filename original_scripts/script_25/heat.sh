#!/bin/bash
#SBATCH --account=brave
#SBATCH --partition=gpu-h100 
#SBATCH --nodes=1
#SBATCH --gres=gpu:h100:1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --job-name=brave
#SBATCH --mem=10GB

module purge
module load mamba cuda/12.2
source /projects/robustmicrob/jlaw/tools/amber24/amber.sh


pmemd.cuda -O -i ../../../script_25/mdt.in -o heat.out -p complex_solvated.prmtop -c min2.rst -r heat.rst -x heat.nc -inf heat.mdinfo -ref min2.rst
