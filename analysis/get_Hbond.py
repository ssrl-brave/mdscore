import mdtraj as md
import numpy as np
import glob

pdb_files = sorted(glob.glob("unique.c0.pdb"))

if not pdb_files:
    raise FileNotFoundError("No PDB files found in the range unique.c0.pdb to unique.c9.pdb.")
def count_hbonds(pdb_file):
    traj = md.load(pdb_file)
    hbonds = md.baker_hubbard(traj, freq=0.2, periodic=False)
    
    # Identify protein and ligand atoms
    top = traj.topology
    protein_atoms = {a.index for a in top.atoms if a.residue.is_protein}
    ligand_atoms = {a.index for a in top.atoms if not a.residue.is_protein}
    
    # Filter H-bonds to only those between protein and ligand
    ligand_hbonds = [hbond for hbond in hbonds if 
                     (hbond[0] in protein_atoms and hbond[2] in ligand_atoms) or 
                     (hbond[0] in ligand_atoms and hbond[2] in protein_atoms)]
    
    return len(ligand_hbonds)

results = []

for pdb in pdb_files:
    hbond_count = count_hbonds(pdb)
    results.append((pdb, hbond_count))
    print(f"{pdb}: {hbond_count} H-bonds")

# Save to a text file
with open("../hbond_counts_boltz.txt", "w") as f:
    for pdb, count in results:
        f.write(f"{pdb} {count}\n")

print("Results saved in hbond_counts.txt")
