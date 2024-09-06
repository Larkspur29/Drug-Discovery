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
import streamlit as st
from padelpy import from_smiles
import pickle
import sklearn


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
    if X is not None and not X.empty:
        X = X['SMILES']
        st.write(X)
        X_np = np.asarray(X)
        a = []
        for i in X_np:
            a.append(i)
        btn1 = st.button('Create Descripter',key=1)
        global res
        global bioactivity_threshold
        if choose == "Corona":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                model = pickle.load(open('corona_model.pkl','rb'))
                res = model.predict(input_x)
                bioactivity_threshold = []
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "CoronaVirus"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')

        elif choose == "Hepatitis A":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                model = pickle.load(open('hepatitisC_model.pkl','rb'))
                res = model.predict(input_x)
                bioactivity_threshold = []
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "Hepatitis C"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')

        elif choose == "Hepatitis B":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                bioactivity_threshold = []
                model = pickle.load(open('hepatitisB_model.pkl','rb'))
                res = model.predict(input_x)
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "Hepatitis B"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')
        elif choose == "Aromatase breast cancer":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                bioactivity_threshold = []
                model = pickle.load(open('aromatase_model.pkl','rb'))
                res = model.predict(input_x)
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "Aromatase"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')


        elif choose == "Dengue":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                model = pickle.load(open('dengue_model.pkl','rb'))
                res = model.predict(input_x)
                bioactivity_threshold = []
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
            df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
            st.markdown("""### Predicted IC50 VALUES""")
            st.write(df)
            os.remove('descripter.csv')
            st.session_state["value"] = "Dengue"
            st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')

        
        elif choose == "Influenza A":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                model = pickle.load(open('influenzaA_model.pkl','rb'))
                res = model.predict(input_x)
                bioactivity_threshold = []
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "Influenza A"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')

        elif choose == "Influenza B":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                model = pickle.load(open('influenzaB_model.pkl','rb'))
                res = model.predict(input_x)
                bioactivity_threshold = []
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "Influenza B"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')

        elif choose == "Alzheimer":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                model = pickle.load(open('alzheimer_model.pkl','rb'))
                res = model.predict(input_x)
                bioactivity_threshold = []
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "Alzheimer"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')


        elif choose == "Chronic kidney":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                model = pickle.load(open('kidney_model.pkl','rb'))
                res = model.predict(input_x)
                bioactivity_threshold = []
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "Chronic Kidney"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')	


        elif choose == "HIV":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                model = pickle.load(open('hiv_model.pkl','rb'))
                res = model.predict(input_x)
                bioactivity_threshold = []
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "HIV"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')		



    #HIV","Leukemia","Lung Cancer	
        elif choose == "LEUKEMIA":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                model = pickle.load(open('leukemia_model.pkl','rb'))
                res = model.predict(input_x)
                bioactivity_threshold = []
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "Leukemia"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')		


    #HIV","Leukemia","Lung Cancer	
        elif choose == "Lung cancer":
            if btn1:
                from_smiles(a,output_csv='descripter.csv',fingerprints=True,descriptors=False)
                input_x = pd.read_csv('descripter.csv')
                input_x = input_x.drop(columns=['Name'],axis=1)
                st.write(input_x)
                model = pickle.load(open('lungs_model.pkl','rb'))
                res = model.predict(input_x)
                bioactivity_threshold = []
                for i in res:
                    if np.log10(i) >= np.log10(10000):
                        bioactivity_threshold.append("inactive")
                    elif np.log10(i) <= np.log10(1000):
                        bioactivity_threshold.append("active")
                    else:
                        bioactivity_threshold.append("intermediate")
                df = pd.concat([X,pd.DataFrame(res),pd.DataFrame({"class":bioactivity_threshold})],axis=1)
                st.markdown("""### Predicted IC50 VALUES""")
                st.write(df)
                os.remove('descripter.csv')
                st.session_state["value"] = "Lung Cancer"
                st.download_button(label="Download data as CSV",data=df.to_csv(),file_name='predictions.csv',mime='text/csv')		


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



    selected_disease = st.selectbox("Select Disease", ["Hepatitis B", "Hepatitis C", "Influenza A",
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
        d_f.to_csv("cleaned_generated_molecules.csv", index=False, columns=["SMILES"])
        X = pd.read_csv("cleaned_generated_molecules.csv")
        file_name = "cleaned_generated_molecules"+".csv"
        st.write("Generated Molecules:")
        sm = pd.read_csv(file_name)
        smiles_list = sm["SMILES"].tolist()
        save_images(smiles_list)
        for smile in smiles_list:
            # Generate 3D structure
            mol = generate_3d_structure(smile)
            if mol is not None:
                img = mol_to_image(mol)
                st.image(img, caption=f'SMILES: {smile}', use_column_width=True)
    X = pd.read_csv("cleaned_generated_molecules.csv")
    screen_molecule(X,selected_disease)
    file_name = "cleaned_generated_molecules"+".csv"
    sm = pd.read_csv(file_name)
    smiles_list = sm["SMILES"].tolist()
    with zipfile.ZipFile('output.zip', 'w') as zipf:
        zipf.write('cleaned_generated_molecules.csv')
        idx=1
        for smile in smiles_list:
            img_path = f'images/image{idx + 1}.png'
        zipf.write(img_path, arcname=os.path.join('images', f'image{idx + 1}.png'))
    with open('output.zip', 'rb') as f:
        st.download_button('Download Zip', f, file_name='output.zip')
# Run the app
if __name__ == "__main__":
    app()
