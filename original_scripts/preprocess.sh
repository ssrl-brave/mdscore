for i in Mpro*
do
        cd $i
        echo $i
        cp Mpro-*protein.pdb receptor_boltz.pdb 
        cp Mpro-*ligand.pdb boltz_lig.pdb
        obabel boltz_lig.pdb -O boltz_lig_h.pdb -h
        grep '^HETATM' boltz_lig_h.pdb>boltz_lig_mod.pdb
        cd ..
done

