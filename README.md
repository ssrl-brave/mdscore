# HTP-MD-ML
Code and instructions for running HTP MD and predict fragment binding. 
As for example, we have here ```Mpro-x0040_iso1``` system under ```mpro_from_boltz``` directory. As for the initial input files, we have boltz generated docked liganded pose ```Mpro-x0040_iso1_model_0_ligand_aligned.pdb``` and receptor ```Mpro-x0040_iso1_model_0_protein.pdb```. 

Step1. Preprocess, parametrize, and submit MD simulations:
---------------------------------------------------------
We first need to add H to the ligands and do a few required formating. We combined all these steps in ```preprocess.sh```. We can simply execute ```bash preprocess.sh``` to generate ```boltz_lig_mod.pdb``` which has a required format for the next steps. 
The next step is to parametrize the ligand and submit MD simulation. This can be performed by executing ```bash generate_inputs_and_submit.sh```. The steps are 
```antechamber -i boltz_lig_mod.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -s 2 -at gaff2 -nc 0
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod
tleap -f ../../script_25/tleap.in
tleap -f ../../script_25/tleap_0.in
parmed -i ../../script_25/hmass.in
bash ../../script_25/chain.sh```
