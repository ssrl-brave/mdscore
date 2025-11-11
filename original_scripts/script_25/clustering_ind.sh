#!/bin/bash
#SBATCH --account=brave
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=shared
#SBATCH --time=12:00:00
#SBATCH --job-name=min
#SBATCH --out=slurm%j.out
#SBATCH --mem-per-cpu=10GB


module purge
module load mamba cuda/12.2
source /projects/robustmicrob/jlaw/tools/amber24/amber.sh
mkdir Clustering_boltz
cd Clustering_boltz


#cluster
cpptraj ../trial_1_boltz/protein.top
cat<<EFO>cpptraj.in

trajin ../trial_1_boltz/protein.nc 1 5000 
trajin ../trial_2_boltz/protein.nc 1 5000 
trajin ../trial_3_boltz/protein.nc 1 5000 

rms first :1-304@CA,CB  out trajrmsd.dat

cluster hieragglo epsilon 2.0 linkage rms :305 nofit sieve 10 summary summary singlerepout representative repout unique repfmt pdb clusterout clusttraj clusterfmt netcdf avgout avg avgfmt pdb out frame_vs_cluster.txt
go
EFO

cpptraj -p ../trial_1_boltz/protein.top -i cpptraj.in
cd ../

cpptraj -i ../../script_25/cpptraj_rmsd_cluster.in

mkdir mmgbsa_0_boltz

cp trial_1_boltz/system_hmass.top mmgbsa_0_boltz

cd mmgbsa_0_boltz

ante-MMPBSA.py -p system_hmass.top -m :1-304 -s :WAT,Na+,Cl- -c com.prmtop -r rec.prmtop -l ligand.prmtop --radii mbondi3
#change the executable path
python /projects/robustmicrob/jlaw/tools/amber24/bin/MMPBSA.py -O -i ../../script_25/mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat -cp com.prmtop -rp rec.prmtop -lp ligand.prmtop -y ../Clustering_boltz/clusttraj.c0

cd ../
