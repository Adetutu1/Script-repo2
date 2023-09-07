# filter compounds based on physico-chemical properites 

from rdkit import Chem

def filter_sdf(input_file, output_file, min_mw, max_mw, min_heavy_atoms, max_heavy_atoms,
               min_aromatic_rings, max_aromatic_rings, min_chiral_centers, max_chiral_centers,
               min_rotatable_bonds, max_rotatable_bonds, min_clogp, max_clogp, min_hbond_acceptors,
               max_hbond_acceptors, min_hbond_donors, max_hbond_donors, min_polar_surface_area,
               max_polar_surface_area, min_num_atoms, max_num_atoms):
    writer = Chem.SDWriter(output_file)

    suppl = Chem.SDMolSupplier(input_file)
    for mol in suppl:
        if mol is not None:
            mw = Chem.Descriptors.MolWt(mol)
            heavy_atoms = Chem.Descriptors.HeavyAtomCount(mol)
            aromatic_rings = Chem.Descriptors.NumAromaticRings(mol)
            chiral_centers = Chem.Descriptors.NumChiralCenters(mol)
            rotatable_bonds = Chem.Descriptors.NumRotatableBonds(mol)
            clogp = Chem.Descriptors.MolLogP(mol)
            hbond_acceptors = Chem.Descriptors.NumHAcceptors(mol)
            hbond_donors = Chem.Descriptors.NumHDonors(mol)
            polar_surface_area = Chem.Descriptors.TPSA(mol)
            num_atoms = mol.GetNumAtoms()

            if (min_mw is None or mw >= min_mw) and (max_mw is None or mw <= max_mw) and \
               (min_heavy_atoms is None or heavy_atoms >= min_heavy_atoms) and \
               (max_heavy_atoms is None or heavy_atoms <= max_heavy_atoms) and \
               (min_aromatic_rings is None or aromatic_rings >= min_aromatic_rings) and \
               (max_aromatic_rings is None or aromatic_rings <= max_aromatic_rings) and \
               (min_chiral_centers is None or chiral_centers >= min_chiral_centers) and \
               (max_chiral_centers is None or chiral_centers <= max_chiral_centers) and \
               (min_rotatable_bonds is None or rotatable_bonds >= min_rotatable_bonds) and \
               (max_rotatable_bonds is None or rotatable_bonds <= max_rotatable_bonds) and \
               (min_clogp is None or clogp >= min_clogp) and \
               (max_clogp is None or clogp <= max_clogp) and \
               (min_hbond_acceptors is None or hbond_acceptors >= min_hbond_acceptors) and \
               (max_hbond_acceptors is None or hbond_acceptors <= max_hbond_acceptors) and \
               (min_hbond_donors is None or hbond_donors >= min_hbond_donors) and \
               (max_hbond_donors is None or hbond_donors <= max_hbond_donors) and \
               (min_polar_surface_area is None or polar_surface_area >= min_polar_surface_area) and \
               (max_polar_surface_area is None or polar_surface_area <= max_polar_surface_area) and \
               (min_num_atoms is None or num_atoms >= min_num_atoms) and \
               (max_num_atoms is None or num_atoms <= max_num_atoms):
                writer.write(mol)

    writer.close()

# Example usage:
filter_sdf('ligprep_7nln_ligands-out.sdf', 'ligprep_7nln_ligands_filter.sdf', min_mw=150, max_mw=500, min_heavy_atoms=, max_heavy_atoms=30,
           min_aromatic_rings=1, max_aromatic_rings=3, min_chiral_centers=0, max_chiral_centers=3,
           min_rotatable_bonds=1, max_rotatable_bonds=6, min_clogp=-1, max_clogp=3, min_hbond_acceptors=1,
           max_hbond_acceptors=8, min_hbond_donors=1, max_hbond_donors=5, min_polar_surface_area=20,
           max_polar_surface_area=130, min_num_atoms=5, max_num_atoms=50)
