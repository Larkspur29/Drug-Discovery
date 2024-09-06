import streamlit as st
import os
from padelpy import from_smiles
import numpy as np
import pickle
import pandas as pd
from PIL import Image

def app():
		file = st.file_uploader("Upload",type=["csv"])
		if file:
			X = pd.read_csv(file)
			X = X['canonical_smiles']
			st.write(X)
			X_np = np.asarray(X)
			a = []
			for i in X_np:
				a.append(i)
			choose = st.selectbox("Choose", ("Hepatitis C","CoronaVirus","Hepatitis B","Influenza A","Influenza B","Aromatase","Alzheimer","Chronic Kidney","Dengue","HIV","Leukemia","Lung Cancer"), key=2)
			btn1 = st.button('Create Descripter',key=1)
			
			global res
			global bioactivity_threshold
			if choose == "CoronaVirus":
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

			elif choose == "Hepatitis C":
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

			elif choose == "Aromatase":
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


			elif choose == "Chronic Kidney":
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

#HIV","Leukemia","Lung Cancer	
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
			elif choose == "Leukemia":
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
			elif choose == "Lung Cancer":
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