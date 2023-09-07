# Filtering compounds using molecular descriptor (physico-chemical properties)

To filter the large compounds dataset based on Lipinski rule of five and other physico-chemical properties parameters is good way to cut down/ filter the large compounds.

To do this, use the 'filter.py' script and adjust the parameters as suitable.

# Filter using Pan Assay Interference (PAINS)

Large compounds can also be filter using PAINS. PAINS are compounds that interferes with biological assays. There are compounds with PAINS pattern and may be eliminated before virtual screening. This would help to eliminate compounds that may interfere with assays if the compounds are taken further for invitro assays.

To do this, use the 'filter_PAINS'.py script.

# Convert sdf files to SMILES

To use the 'filter_PAINS' script the single compound file must be converted to SMILES in CSV format with two columns; one column represents the ZINC ID while the other column represents each compounds SMILES.

Use the 'convert_sdf_to_smile.py' script to achieve this.