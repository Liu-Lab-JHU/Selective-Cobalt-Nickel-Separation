#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Morgan Fingerprints generation
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

# Load the CSV file containing SMILES
df = pd.read_csv('file.csv')         # Replace with the path to your CSV file
smiles_list = df['SMILES'].tolist()  # Replace 'SMILES' with the actual column name

fingerprint_data = []

# Function to convert SMILES to RDKit molecule and then to Morgan Fingerprint
def smiles_to_morgan(smiles, radius=2, nBits=2048):
    """
    :param smiles: str, SMILES string of the molecule
    :param radius: int, radius of the Morgan fingerprint. Default is 2
    :param nBits: int, number of bits in the Morgan fingerprint. Default is 2048
    :return: Morgan fingerprint of the molecule
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    else:
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)

# Apply the function to each SMILES in the list
for smi in smiles_list:
    morgan_fps = smiles_to_morgan(smi)
    bit_list = list(morgan_fps.ToBitString())
    fingerprint_data.append(bit_list)

columns = list(range(2048))

# Optional: Save the fingerprints to a new CSV file
df = pd.DataFrame(fingerprint_data, columns=columns)
df.to_csv('morganFP.csv', index=False)