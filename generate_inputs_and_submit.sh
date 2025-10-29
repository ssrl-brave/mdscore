antechamber -i boltz_lig_mod.pdb -fi pdb -o ligand.mol2 -fo mol2 -c bcc -s 2 -at gaff2 -nc 0
parmchk2 -i ligand.mol2 -f mol2 -o ligand.frcmod
tleap -f ../../script_25/tleap.in
tleap -f ../../script_25/tleap_0.in
parmed -i ../../script_25/hmass.in
bash ../../script_25/chain.sh
