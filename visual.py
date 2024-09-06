import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
from stmol import showmol

# Function to read SMILES from CSV file
def read_smiles_from_csv(file):
    df = pd.read_csv(file)
    return df['SMILES'].tolist()

# Streamlit app
def app():
    st.title("Molecular Visualization")
    
    # File uploader for CSV
    uploaded_file = "cleaned_generated_molecules.csv"
    
    if uploaded_file is not None:
        smiles_list = read_smiles_from_csv(uploaded_file)
        
        # Dropdown for selecting SMILES
        selected_smiles = st.selectbox("Select SMILES", smiles_list)
        
        # Display molecular model
        if selected_smiles:
            mol = Chem.MolFromSmiles(selected_smiles)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
            mblock = Chem.MolToMolBlock(mol)
            width = 700
            height = 700
            cartoon_radius = 0.2
            stick_radius = 0.2
            style = "stick"
            view = py3Dmol.view(width=width, height=height)
            view.addModel(mblock, 'mol')
            view.setStyle({style:{}})
            view.spin()
            showmol(view, height=height, width=width)
