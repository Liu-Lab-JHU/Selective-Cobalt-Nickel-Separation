#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Feature engineering
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

"""
Step 1: generation of the ECFPs for the bioacids.
"""
# Load the CSV file containing SMILES
df = pd.read_csv('raw_data.csv')         # Replace with the path to your CSV file
smiles_list = df['SMILES'].tolist()  # Replace 'SMILES' with the actual column name

fingerprint_data = []

# Function to convert SMILES to RDKit molecule and then to ECFP
def smiles_to_morgan(smiles, radius=2, nBits=2048):
    """
    :param smiles: str, SMILES string of the molecule
    :param radius: int, radius of the ECFP. Default is 2
    :param nBits: int, number of bits in the ECFP. Default is 2048
    :return: ECFP of the molecule
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
df.to_csv('dataset_w_ECFPs.csv', index=False)

"""
Step 2: generation of the WECFPs for the bioacids.
To take the concentration of the bioacids into account, we used normalized concentration as the weight 
of the ECFPs to generate the weighted ECFPs (denoted by WECFPs).
"""
feature_columns = [col for col in df.columns if col not in ['acid','SMILES','concentration','potential','selectivity']]
weight_column = 'concentration'
normalized_weight=60

# Multiply each feature by the weight column
df[feature_columns] = df[feature_columns].mul(df[weight_column], axis=0).div(normalized_weight)

# Optionally, you can save the resulting dataset
df.to_csv('dataset_w_MorganFPs_weighted.csv', index=False)