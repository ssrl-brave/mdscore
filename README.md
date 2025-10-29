# HTP-MD-ML
Code and instructions for running HTP MD and predict fragment binding. 
As for example, we have here ```Mpro-x0040_iso1``` system under ```mpro_from_boltz``` directory. As for the initial input files, we have boltz generated docked liganded pose ```Mpro-x0040_iso1_model_0_ligand_aligned.pdb``` and receptor ```Mpro-x0040_iso1_model_0_protein.pdb```. 

Step1. Preprocess, parametrize, and submit MD simulations:
---------------------------------------------------------
We first need to add H to the ligands and do a few required formating. We combined all these steps in ```preprocess.sh```. We can simply execute ```bash preprocess.sh``` to generate ```boltz_lig_mod.pdb``` which has a required format for the next steps. 
The next step is to parametrize the ligand and submit MD simulation. This can be performed by executing ```bash generate_inputs_and_submit.sh```. The steps are 
```
antechamber -i boltz_lig_mod.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -s 2 -at gaff2 -nc 0
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod
tleap -f ../../script_25/tleap.in
tleap -f ../../script_25/tleap_0.in
parmed -i ../../script_25/hmass.in
bash ../../script_25/chain.sh
```
The first three lines are getting the amber parameter for the ligand. The net charge (nc) is 0 here but should be changed depending on the system. Once we have the parameter we execute ```tleap_0.in``` to get the solvated protein-ligand complex parameter and coordinates which are the inputs for running MD simulation. ```hmass.in``` file is used to perform hydrogen mass repartition that allowes us to use higher timestep in MD. Finally by executing ```chain.sh``` we submit MD simulation job where each steps are depended on the next step. To elaborate, the ```chain.sh``` has
```
#!/bin/bash

mkdir trial_1_boltz trial_2_boltz trial_3_boltz

script=../../script_25/min.sh

PID0=$(sbatch $script|awk '{print $NF}')
echo $PID0

for i in trial_1_boltz trial_2_boltz trial_3_boltz 
do
cd $i

script=../../../script_25/heat.sh

PID_1=$(sbatch --dependency=afterany:$PID0  $script| awk '{print $NF}')
echo $PID_1

script=../../../script_25/npt_cpu.sh

PID_2=$(sbatch --dependency=afterany:$PID_1  $script| awk '{print $NF}')
echo $PID_2

script=../../../script_25/md.sh

PID_3=$(sbatch --dependency=afterany:$PID_2  $script| awk '{print $NF}')
echo $PID_3

cd ..
done
```

We first create 3 directories for running the simulation in triplicate. We then minimize the generated initial coordinates and then perform heating, cpu based equilbration and then finally gpu based equilibration and production. All input files are in ```script_25``` directory. 

Step2. Post processing of the trajectory and featurization:
----------------------------------------------------------

Once the MD simulations are done, we perform a series of analysis---
a. ```cpptraj -i ../../script_25/cpptraj_strip.in``` run it under ```trial_*_boltz``` directory to generate the dry docked complex that is used as the reference structure. 

b. ```cpptraj -i ../../script_25/cpptraj_strip_2.in``` run it under ```trial_*_boltz``` directory to strip off waters and ions from the generated trajectory and parameter. This gives us dry trajectory ```protein.nc``` and parameter ```protein.top```.

c. ```cpptraj -i ../../script_25/cpptraj_rmsd.in``` run it under ```trial_*_boltz``` directory to calcualte ligand and protein rmsd of the trajectory taking the docked structure as reference. 
d. ```cpptraj -i ../../script_25/cpptraj_rmsf.in``` run it under ```trial_*_boltz``` directory to calcualte ligand rmsf over the trajectory. 
e. ```sbatch -i ../../script_25/clustering_ind.sh``` run it under base system directory to cluster the ligand poses based on their binding similarity. Once, clustering is done, this script also calculate RMSD of the top 5 cluster centers and perform MMGBSA energy calculation on the top cluster. 
f. ```cpptraj -i ../../script_25/get_water.in``` run it under ```trial_*_boltz``` directory this uses solvated trajectory to calculate number of water molecule in the active site during the simulation.
g. ```python ../../analysis/get_Hbond.py``` run it under ```clustering_boltz``` directory to estimate number of H bonds between ligand and receptor in top representative structure from MD.
h. ```python ../../analysis/get_distance.py``` run it under ```clustering_boltz``` directory to estimate minumum average distance between the center of mass of ligand and key anchoring residues of receptor in the top representative structure from MD.
i. ```python ../../analysis/get_distance_boltz.py``` run it under base system directory to estimate the same distances for the docked structure. 
j. ```python ../../analysis/get_energies.py``` run it outside the base system directies to grab the MMGBSA energy values. 
k. ```python ../../analysis/get_affinity.py``` run in the base system directories to grab the affinity predicted by boltz docking. 

and finally,
l. ```python ../../analysis/get_all_feature_for_ML_boltz.py``` run it outside the base system directies to combine all these values and write a csv file. 

Step3. ML models:
----------------------------------------------------------
