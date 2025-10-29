import gemmi

# Load the CIF file
cif_file = "Mpro-x0040_iso1_model_0.cif"
structure = gemmi.read_structure(cif_file)

# Save as PDB
pdb_file = "boltz_modeled.pdb"
structure.write_pdb(pdb_file)

print(f"Converted {cif_file} to {pdb_file} including ligands.")
