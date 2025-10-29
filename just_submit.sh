#!/bin/bash
#SBATCH --account=brave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared
#SBATCH --time=2:00:00
#SBATCH --job-name=water
#SBATCH --out=slurm%j.out
#SBATCH --mem-per-cpu=1GB

cpptraj -i ../../script_25/get_water.in
