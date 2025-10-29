from rdkit import Chem
from rdkit.Chem import AllChem
sml=open('out.smi').readline().strip()
mol_from_smiles = Chem.MolFromSmiles(sml)
#mol_with_h = Chem.AddHs(mol_from_smiles)

template_mol = Chem.MolFromPDBFile('Mpro-x0040_iso1_model_0_ligand.pdb', removeHs=False)

aligned_mol = AllChem.AssignBondOrdersFromTemplate( mol_from_smiles,template_mol)

#aligned_mol_h=Chem.AddHs(aligned_mol)
#Chem.MolToPDBFile(aligned_mol_h, "aligned_with_H.pdb")
Chem.MolToPDBFile(aligned_mol, "Mpro-x0040_iso1_model_0_ligand_aligned.pdb")

