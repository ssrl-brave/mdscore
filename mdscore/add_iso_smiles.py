from rdkit import Chem
from rdkit.Chem import EnumerateStereoisomers

def add_iso_smiles(df_):
    df = df_.copy()
    assert "Smiles" in df
    all_iso_smiles = []
    for i_smi, smi in enumerate(df.Smiles):
        mol = Chem.MolFromSmiles(smi)
        options = EnumerateStereoisomers.StereoEnumerationOptions(unique=True, tryEmbedding=True)
        isomers = tuple(EnumerateStereoisomers.EnumerateStereoisomers(
            mol,
            options=options)
        )

        # prepare molecules further:
        isomers = [Chem.AddHs(mol) for mol in isomers]
        # end prepare molecules further

        print(f"{i_smi + 1}/{len(df)} found {len(isomers)} smiles strings for fraser smiles {smi}")
        iso_smiles = []
        for i_mol, mol in enumerate(isomers):
            mol = Chem.AddHs(mol)

            iso_smile = Chem.MolToSmiles(mol, isomericSmiles=True)
            iso_smiles.append(iso_smile)
            #mol_id = f"base{i_smi}-iso{i_mol}"
            #sdf = f"{outdir}/{mol_id}.sdf"
            #mol.SetProp("ID", mol_id)
            #with Chem.SDWriter(sdf) as writer:
            #    writer.write(mol)

            #sdf_record = f"{sdf}:{iso_smile}"
            #sdf_records.append(sdf_record)
        all_iso_smiles.append(";".join(iso_smiles))

    #df["sdf_records"] = all_sdf_records
    df["iso_smiles"] = all_iso_smiles
    df = df.assign(iso_smiles=df.iso_smiles.str.split(';')).explode('iso_smiles').reset_index(drop=True)
    #df[["rdkit_sdf", "rdkit_isomer_smiles"]] = df.sdf_records.str.split(":", n=1, expand=True)
    #df.to_csv(f"{csv_f}", sep="\t", index=False)
    return df
