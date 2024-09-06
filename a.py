from transformers import RobertaTokenizer, RobertaForCausalLM
import streamlit as st
import csv
import os
import shutil
import time
from padelpy import padeldescriptor
import numpy as np
import pandas as pd
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from sklearn.feature_selection import VarianceThreshold
import zipfile



# Function to generate 3D structure from SMILES
def generate_3d_structure(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol)
        return mol
    else:
        return None

# Function to convert RDKit molecule to image
def mol_to_image(mol):
    # Convert RDKit molecule to image
    img = Draw.MolToImage(mol, size=(400, 400))
    return img

def generate_similar_mol(input_sequences, num_samples=10):
    tokenizer = RobertaTokenizer.from_pretrained("seyonec/PubChem10M_SMILES_BPE_450k", padding_side="left")
    model = RobertaForCausalLM.from_pretrained("seyonec/PubChem10M_SMILES_BPE_450k", is_decoder=True)
    
    sim_mol = set()
    created_data = []
    
    for input_smiles in input_sequences:
        try:
            input_ids = tokenizer.encode(input_smiles, return_tensors="pt")
            outputs = model.generate(input_ids, max_length=60, num_return_sequences=num_samples, temperature=0.7, do_sample=True)
            
            for output_ids in outputs:
                decoded_sequence = tokenizer.decode(output_ids, skip_special_tokens=True)
                mol = Chem.MolFromSmiles(decoded_sequence)
                if mol is not None:
                    mol_smiles = Chem.MolToSmiles(mol)
                    
                    if mol_smiles not in sim_mol:
                        AllChem.Compute2DCoords(mol)
                        sim_mol.add(mol_smiles)
        except Exception as e:
            st.error(f"Failed to process molecule: {input_smiles}, Error: {str(e)}")
    return pd.DataFrame(sim_mol, columns=["SMILES"])
# Function to generate similar molecules from input SMILES

# Function to read SMILES from CSV file
def read_smiles_from_csv(file_path):
    smiles_list = []
    with open(file_path, 'r', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip header
        for row in reader:
            if row and len(row) > 0:
                smiles_list.append(row[0])
    return smiles_list

# Function to filter molecules based on bioactivity threshold
def screen_molecule(X, choose):
    global bioactivity_threshold
    d = pd.read_csv("cleaned_generated_molecules.csv")
    d2 = pd.concat([d["SMILES"]], axis=1)
    d2.to_csv('molecule.smi', sep='\t', index=False, header=False)
    padeldescriptor(mol_dir='molecule.smi',
                    d_file='Substructure.csv',
                    detectaromaticity=True,
                    standardizenitro=True,
                    standardizetautomers=True,
                    threads=2,
                    removesalt=True,
                    log=True,
                    fingerprints=True)
    descriptors = pd.read_csv('Substructure.csv')
    descriptors.drop(index=0, columns=['Name'], axis=1, inplace=True)
    bioactivity_threshold = []
    for i in range(len(descriptors)):
        s = descriptors.iloc[i, :].sum()
        log_value = np.log10(s)
        if log_value >= np.log10(10000):
            bioactivity_threshold.append("inactive")
        elif log_value <= np.log10(1000):
            bioactivity_threshold.append("active")
        else:
            bioactivity_threshold.append("intermediate")
    df = pd.concat([X, pd.DataFrame({"class": bioactivity_threshold})], axis=1)
    df_active = df[df['class'] == 'active']
    s = choose+".csv"
    df_active.to_csv(s, index=False)
    return df_active


def save_images(smiles_list):
    os.makedirs('images', exist_ok=True)
    for idx, smile in enumerate(smiles_list):
        mol = generate_3d_structure(smile)
        if mol is not None:
            img = mol_to_image(mol)
            img_path = f'images/image{idx + 1}.png'  # Use generic names like "image1", "image2", etc.
            img.save(img_path)


# Streamlit app
# Streamlit app
def app():
    st.title("Molecule Generator")



    selected_disease = st.selectbox("Select Disease", ["Hepatitis B", "Hepatitis A", "Influenza A",
                                                    "Influenza B", "Aromatase breast cancer",
                                                    "Alzheimer", "Dengue", "Chronic kidney",
                                                    "Lung cancer", "Corona", "HIV", "LEUKEMIA"])
    sample_file_mapping = {
        "Hepatitis B": "hepatitis_b.csv",
        "Hepatitis C": "hepatitis_c.csv",
        "Influenza A": "influenza_a.csv",
        "Influenza B": "influenza_b.csv",
        "Aromatase breast cancer": "aromatase_breast_cancer.csv",
        "Alzheimer": "alzheimer.csv",
        "Dengue": "dengue.csv",
        "Chronic kidney": "chronic_kidney.csv",
        "Lung cancer": "lung_cancer.csv",
        "Corona": "corona.csv",
        "HIV": "hiv.csv",
        "LEUKEMIA": "leukemia.csv"
    }

    if st.button("Generate",key="hello"):
        selected_file = sample_file_mapping[selected_disease]
        samplemol = read_smiles_from_csv(selected_file)
        newmol = generate_similar_mol(samplemol)
        newmol.to_csv("generated_molecules.csv")
        d = pd.read_csv("generated_molecules.csv")
        special_characters = "@#!$%^&*"
        d['SMILES'] = d['SMILES'].replace(f'[{special_characters}]', '', regex=True)
        d_f = d[d['SMILES'].str.len() > 6]
        d_f.to_csv("cleaned_generated_molecules.csv", index=False)
        X = pd.read_csv("cleaned_generated_molecules.csv")
        df_active = screen_molecule(X, selected_disease)
        file_name = selected_disease+".csv"
        st.write("Generated Molecules:")
        st.write(df_active)
        sm = pd.read_csv(file_name)
        smiles_list = sm["SMILES"].tolist()
        save_images(smiles_list)
        for smile in smiles_list:
            # Generate 3D structure
            mol = generate_3d_structure(smile)
            if mol is not None:
                img = mol_to_image(mol)
                st.image(img, caption=f'SMILES: {smile}', use_column_width=True)
    file_name = selected_disease+".csv"
    sm = pd.read_csv(file_name)
    smiles_list = sm["SMILES"].tolist()
    with zipfile.ZipFile('output.zip', 'w') as zipf:
        zipf.write('cleaned_generated_molecules.csv')
        for idx, smile in smiles_list:
            img_path = f'images/image{idx + 1}.png'
        zipf.write(img_path, arcname=os.path.join('images', f'image{idx + 1}.png'))
    with open('output.zip', 'rb') as f:
            st.download_button('Download Zip', f, file_name='output.zip')

# Run the app
if __name__ == "__main__":
    app()
