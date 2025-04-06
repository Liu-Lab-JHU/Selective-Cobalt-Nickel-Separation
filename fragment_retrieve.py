#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Retrieve and draw the fragments for each important bit (or specific bits)
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import rdMolDraw2D

df = pd.read_csv('raw_data.csv')  # Replace with the path to your CSV file that contains SMILES of molecules
smiles_list = df['SMILES'].tolist() # Replace 'SMILES' with the actual column name
molecules = []
fingerprint_data = []
nBits=2048  # Replace to actual number of bits

# Function to convert SMILES to RDKit molecule and then to ECFP
for smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)  # Convert SMILES to RDKit molecule
    if mol:
        # Generate Morgan Fingerprint
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=nBits)
        # Convert to bit list
        bit_list = list(fp.ToBitString())
        fingerprint_data.append(bit_list)
        molecules.append(mol)  # Store the molecule for later use
    else:
        print(f"Invalid SMILES: {smiles}")
        fingerprint_data.append([None] * nBits)


df_fragments=pd.read_csv('feature_importance.csv')  # Replace with the path to your CSV file
important_bits=[]

for feature in df_fragments['feature']:
    important_bits.append(feature)

# important_bits = [1141, 216, 1694, 1564]  # Replace with important bits
# Generate Morgan FP with bitInfo
info = {}


for mol in molecules:  # Use stored molecules
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048, bitInfo=info)

    # Check fragments for each important bit
    for bit in important_bits:
        if bit in info:
            print(f"\nFragment for bit {bit}:")
            atom_idx, radius = info[bit][0]  # Take the first associated atom environment
            env = Chem.FindAtomEnvironmentOfRadiusN(mol, 2, atom_idx)  # Get atom environment
            amap = {}
            submol = Chem.PathToSubmol(mol, env, atomMap=amap)  # Substructure for this bit

            # Draw the substructure fragment and save it as an SVG file
            drawer = rdMolDraw2D.MolDraw2DSVG(400, 400)
            rdMolDraw2D.PrepareAndDrawMolecule(drawer, submol)
            drawer.FinishDrawing()
            with open(f"fragment_bit_{bit}.svg", "w") as f:
                f.write(drawer.GetDrawingText())
            print(f"Fragment for bit {bit} saved as 'fragment_bit_{bit}.svg'")
        else:
            print(f"Bit {bit} not found in this molecule.")