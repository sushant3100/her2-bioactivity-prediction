# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:49:29 2021

@author: Sushant
"""
import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle
import glob
from padelpy import padeldescriptor

xml_files = glob.glob('*.xml')

#Creating a list of present files
FP_list = ['PubChem']

#Creating Data Dictionary
fp = dict(zip(FP_list, xml_files))
# Molecular descriptor calculator
def padel_desc():
    """This function will calculate Molecular Fingerprints 
        based on PaDEL descriptors"""
    fingerprints = 'PubChem'
    
    fingerprint_output_file=''.join([fingerprints,'.csv'])
    fingerprint_descriptortypes = fp[fingerprints]
    
    padeldescriptor(mol_dir='molecule.smi',
                     d_file=fingerprint_output_file,
                     descriptortypes=fingerprint_descriptortypes,
                     detectaromaticity=True,
                     standardizenitro=True,
                     standardizetautomers=True,
                     threads=2,
                     removesalt=True,
                     log=True,
                     fingerprints=True)
    
# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href
    
 # Model building
def build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('HER2_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = pd.Series(load_data[1], name='molecule_name')
    df = pd.concat([molecule_name, prediction_output], axis=1)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)  
    
#logo image
image = Image.open('logo.jpg')
 
st.image(image, use_column_width=True)

#Page title
st.markdown("""
            # Bioactivity Prediction App for HER2 inhibitors
            
            This app allows you to predict the bioactivity of input compounds towards inhibiting the 'HER2' protien which is found on the cancer cells in breast cancer.
            
            **Note**: If pIC50 value is greater then 6 the bioactivity for that compound is **'active'**.If it is less than 5, bioactivity is **'inactive**.
            Bioactivity is stated as intermediate if the value lies between 5 and 6.
            """) 
#sidebar
with st.sidebar.header('1.Upload your CSV data'):
    uploaded_file=st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
                        [Example input file](https://raw.githubusercontent.com/dataprofessor/bioactivity-prediction-app/main/example_acetylcholinesterase.txt)
""")
    
if st.sidebar.button('Predict'):
    load_data = pd.read_table(uploaded_file, sep=' ', header=None)
    load_data.to_csv('molecule.smi',sep='\t', header=False, index = False)
    
    st.header('**Original input data**')
    st.write(load_data)
    
    with st.spinner("Calculating descriptors please wait....."):
        padel_desc()
        
    # Read in calculated descriptors and display the dataframe
    st.header('**Calculated molecular descriptors**')
    desc = pd.read_csv('PubChem.csv')
    st.write(desc)
    st.write(desc.shape)
    
    # Read descriptor list used in previously built model
    st.header('**Subset of descriptors from previously built models**')
    Xlist = list(pd.read_csv('descriptor_list.csv').columns)
    desc_subset = desc[Xlist]
    st.write(desc_subset)
    st.write(desc_subset.shape)
    
     # Apply trained model to make prediction on query compounds
    build_model(desc_subset)
    
else:
    st.info('Upload input data in the sidebar to start!')