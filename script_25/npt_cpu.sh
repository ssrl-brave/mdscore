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


pmemd -O -i ../../../script_25/npt1_cpu.in -o npt1_cpu.out -p complex_solvated.prmtop -c heat.rst -r npt_cpu.rst -inf npt_cpu.mdinfo -ref heat.rst



