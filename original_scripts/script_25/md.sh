#!/bin/bash
#SBATCH --account=brave
#SBATCH --partition=gpu-h100 
#SBATCH --nodes=1
#SBATCH --gres=gpu:h100:1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=1
#SBATCH --time=36:00:00
#SBATCH --job-name=brave
#SBATCH --mem=10GB

module purge
module load mamba cuda/12.2
source /projects/robustmicrob/jlaw/tools/amber24/amber.sh


pmemd.cuda -O -i ../../../script_25/npt1.in -o npt1.out -p complex_solvated.prmtop -c npt_cpu.rst -r npt1.rst -x npt1.nc -inf npt1.mdinfo -ref npt_cpu.rst
pmemd.cuda -O -i ../../../script_25/npt2.in -o npt2.out -p complex_solvated.prmtop -c npt1.rst -r npt2.rst -x npt2.nc -inf npt2.mdinfo -ref npt1.rst
pmemd.cuda -O -i ../../../script_25/md.in -o md.out -p system_hmass.top -c npt2.rst -r md.rst -x md.nc -inf md.mdinfo -ref npt2.rst
