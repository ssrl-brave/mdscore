import mdtraj as md
import numpy as np
import glob

# Define residue numbers (adjusting for 0-based indexing in MDTraj)
#residue_numbers = [19,21,46,47,50,152,154]  # PDB numbering
residue_numbers = [164,141,165,171,139,162,144,40,48]

# Find all available PDB files (from unique.c0.pdb to unique.c9.pdb)
pdb_files = sorted(glob.glob("unique.c0.pdb"))

if not pdb_files:
    raise FileNotFoundError("No PDB files found in the range unique.c0.pdb to unique.c9.pdb.")

# Open output file to store results
with open("../avg_distances_boltz.txt", "w") as out_file:
    #out_file.write("PDB File\tAverage Distance (nm)\n")

    for pdb_file in pdb_files:
        try:
            # Load PDB
            pdb = md.load(pdb_file)

            # Get residue indices
            residue_indices = [res.index for res in pdb.topology.residues if res.resSeq in residue_numbers]

            # Select ligand atoms (assuming residue name is UNL)
            ligand_atoms = pdb.topology.select("resname LIG")
            if len(ligand_atoms) == 0:
                print(f"Skipping {pdb_file}: No ligand (UNL) found.")
                continue

            # Compute ligand center of mass (COM)
            ligand_com = md.compute_center_of_mass(pdb.atom_slice(ligand_atoms))

            # Compute closest distances for each residue
            closest_distances = []
            for res_idx in residue_indices:
                residue_atoms = [a.index for a in pdb.topology.atoms if a.residue.index == res_idx]
                res_positions = pdb.xyz[0, residue_atoms, :]  # Extract atom positions

                # Compute Euclidean distances to ligand COM
                distances = np.linalg.norm(res_positions - ligand_com[0], axis=1)

                # Store the minimum distance for this residue
                closest_distances.append(np.min(distances))

            # Compute average closest distance
            avg_distance = np.mean(closest_distances)

            # Print and save results
            print(f"{pdb_file}: Average distance = {avg_distance:.3f} nm")
            out_file.write(f"{pdb_file}\t{avg_distance:.3f}\n")

        except Exception as e:
            print(f"Error processing {pdb_file}: {e}")

