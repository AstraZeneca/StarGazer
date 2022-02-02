# StarGazer

![Maturity level-Prototype](https://img.shields.io/badge/Maturity%20Level-Prototype-red)

StarGazer is a tool designed for rapidly assessing drug repositioning opportunities. It combines multi-source, multi-omics data with a novel target prioritization scoring system in an interactive Python-based Streamlit dashboard. StarGazer displays target prioritization scores for genes associated with 1844 phenotypic traits, from type II diabetes and chronic hepatitis B, to hypertension and IgE levels. 

Target prioritization is essential for drug discovery and repositioning. Applying computational methods to analyze and process multi-omics data in order to find new drug targets is a practical approach for achieving this. Despite an increasing number of methods for generating datasets such as genomics, phenomics, and proteomics, attempts to integrate and mine such datasets remain limited in scope. We believe that integrating different data sources using a singular numerical scoring system could help to bridge these different omics layers and facilitate rapid drug target prioritization for studies in drug discovery, development or repositioning. 

 

# How to run StarGazer
  + This guide is for Windows. Please cater these instructions to your OS
  + Python 3.9.7 was used for building StarGazer

1. Download this repository and unzip the folder to your desired location
  + Do not move any of the individual items from this folder
  + (Optional) Create a virtual environment at this stage

2. Open the command line and change the working directory to this folder. This may require using a command that looks like this:

~~~
cd C:\Users\personal\Downloads\StarGazer
~~~

3. Install package dependencies using requirements.txt (may take a few minutes). This may require using a command that looks like this:

~~~
pip install -r requirements.txt
~~~
  + If there are any problem in installing packages, please try directly using pip install, e.g., "pip install streamlit"

4. Run Streamlit on the StarGazer Python script using the following line of code:

~~~
streamlit run StarGazer.py
~~~


