# To convert sdf compounds in a single sdf file into smiles as a CSV file

import pandas as pd
from rdkit import Chem

def convert_sdf_to_smiles_csv(input_file, output_file):
    supplier = Chem.SDMolSupplier(input_file)
    data = []
    for mol in supplier:
        if mol is not None:
            cid = mol.GetProp("_Name")  # Assuming compound ID is stored as a property named "_Name"
            smiles = Chem.MolToSmiles(mol)
            data.append({'CompoundID': cid, 'SMILES': smiles})
    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)

# Provide the paths to your input SDF file and desired output CSV file
input_sdf_file = '/path/to/input/file/'
output_csv_file = '/path/to/output/file'

convert_sdf_to_smiles_csv(input_sdf_file, output_csv_file)
