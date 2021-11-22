from logging import NullHandler
from numpy.__config__ import show
from pkg_resources import yield_lines
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from streamlit_echarts import st_echarts
import requests
import json
from pyvis import network as net
from stvis import pv_static
import io
import collections
from sklearn import preprocessing
import base64
from io import BytesIO
import os
from PIL import Image
import webbrowser

# download function
def get_table_download_link(df):
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    csv = df.to_csv(index=True)
    b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}">Download data</a>'
    return href

# function to generate a data frame of gene symbol and openTargets association score
def opentargets_gene_score(disease_name):

    # Set base URL of GraphQL API endpoint
    base_url = "https://api.platform.opentargets.org/api/v4/graphql"
    
    # query disease id via GraphQL API
    query_string1 = """
        query searchDiseaseID($diseaseName: String!, $entityNames: [String!]) {
          search(queryString: $diseaseName, entityNames: $entityNames ) {
            total
          hits{
            id
            name
          }
          }
        }
    """
    
    query_string2 = """
        query associatedTargets($diseaseID: String!) {
          disease(efoId: $diseaseID) {
            id
            name
            associatedTargets(page: { index: 0, size: 300 })  {
              count
              rows {
                target {
                  approvedSymbol
                }
                score
              }
            }
          }
        }
    """
    # Set variables object of arguments to be passed to endpoint
    variables = {"diseaseName": disease_name, "entityNames": ["disease"]}
    
    # Perform POST request and check status code of response
    r = requests.post(base_url, json={"query": query_string1, "variables": variables})
    
    try:
        df = pd.json_normalize(r.json()["data"]["search"]["hits"])
        disease_id = df.loc[df["name"].str.lower() == disease_name.lower(), "id"].values[0]
        variables = {"diseaseID": disease_id}

        r = requests.post(base_url, json={"query": query_string2, "variables": variables})
        gene_scoreDF = pd.json_normalize(r.json()["data"]["disease"]["associatedTargets"]["rows"])
        gene_scoreDF = gene_scoreDF.rename({
            "score": "opentargets_associations",
            "target.approvedSymbol": "gene_symbol"
        }, axis = 1)
    except:
        gene_scoreDF = []
            
    return gene_scoreDF

def proteins_interaction(input_protein):
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "network"
    request_url = "/".join([string_api_url, output_format, method])
    params = {
    "identifiers" : "%0d".join(input_protein), # your protein in a list
    "species" : 9606, # species NCBI identifier 
    "caller_identity" : "stargazer" # your app name
    }
    response = requests.post(request_url, data=params)
    r = response.content
    rawData = pd.read_csv(io.StringIO(r.decode('utf-8')), sep = "\t")
    rawData_drop = rawData.drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
    return rawData_drop

def go_enrichment(input_gene):
    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "enrichment"
    request_url = "/".join([string_api_url, output_format, method])
    params = {
        "identifiers" : "%0d".join(input_gene), # your protein
        "species" : 9606, # species NCBI identifier 
        "caller_identity" : "stargazer" # your app name
    }
    response = requests.post(request_url, data=params)
    r = response.content
    rawData = pd.read_csv(io.StringIO(r.decode('utf-8')), sep = "\t")
    return rawData

# initialising app
path = os.getcwd()
    
## main page set up
st.set_page_config(layout="wide", page_title="StarGazer")

# import dataset
df = pd.read_csv(path + "/assets/phewas-catalog.csv")

# fill na with "Unknown"
df['gene_name'] = df['gene_name'].fillna("UNKNOWN")
df_selected = df[["gene_name", "snp", "phewas phenotype", "p-value", "odds-ratio", "gwas-associations"]]

# Adding COVID-19 module
full_url = "https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0100096/associations?projection=associationByEfoTrait"
r = requests.get(full_url, json={})
if r.status_code == 200:
    COVID_df = pd.json_normalize(r.json()["_embedded"]["associations"], ["snps", ["genomicContexts"]], ["orPerCopyNum", "pvalue"], errors = "ignore")[["gene.geneName", "pvalue", "_links.snp.href", "orPerCopyNum"]].drop_duplicates(keep = "first")
    COVID_df["_links.snp.href"] = COVID_df["_links.snp.href"].str.strip("{?projection}").str.split("Polymorphisms/").str[1]
    COVID_df = COVID_df.loc[COVID_df["orPerCopyNum"].notna(), :].reset_index(drop = True)
    COVID_df["orPerCopyNum"] = COVID_df["orPerCopyNum"].astype(float)
    COVID_df = COVID_df.rename({
        "gene.geneName": "gene_name",
        "pvalue": "p-value",
        "_links.snp.href": "snp",
        "orPerCopyNum": "odds-ratio"
    }, axis = 1)
    COVID_df[["phewas phenotype", "gwas-associations"]] = "COVID-19"
df_selected = df_selected.append(COVID_df).reset_index(drop = True)

# extract data from Pharos
query_string = """
query AllTargets {
    targets(
        filter: { 
          facets: [{
              facet: "Target Development Level",
                values: ["Tclin", "Tchem", "Tbio", "Tdark"]
            }]
        }
    ) {
        targets (top : 100000) {
            sym
            tdl
        }
    }
}
"""

r = requests.post("https://pharos-api.ncats.io/graphql", json={"query": query_string})
if r.status_code == 200:
    df_druggable = pd.DataFrame(r.json()["data"]["targets"]["targets"]).drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)


# loading logo
col1, col2 = st.columns((4, 1))
path = os.getcwd()
logo = Image.open(path + "/assets/logo.png")
col2.image(logo, output_format = "PNG", width = 200)

# tile of the dashboard and description
st.markdown("***")
st.title("StarGazer: Multi-omics evidence-based drug target prioritization")

select = st.sidebar.selectbox('Search by', ["--", 'Gene', 'Variant', 'PheWAS', 'GWAS', 'GWAS_PheWAS Union', 'GWAS_PheWAS Intersection', "Protein-protein Interaction", 'Disease Target Prioritization'], key='1')



if select == "Gene":
    st.markdown("This dashboard shows the associated phenotypes of your genes of interest.")
    # sidebar -- gene & variant select boxs
    gene = sorted(df_selected["gene_name"].unique().tolist())
    select_gene = st.sidebar.selectbox('Gene', gene, key='2')
    variant = sorted(df_selected[df_selected["gene_name"] == select_gene]["snp"].unique().tolist())
    select_variant = st.sidebar.selectbox('Variant', variant, key='3')

    # subset the data frame
    df_variant = df_selected[df_selected["snp"] == select_variant]

    # sidebar --  p-value slider
    df_variant.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
    select_p = st.sidebar.text_input(label = "P-value", help = "Defaults to p = 0.05. Accepts scientific notation, e.g., 5E-4, 3e-9", value = "0.05")
    try:
        if (float(select_p) <= 1) & (float(select_p) > 0):
            select_p = float(select_p)
        else:
            select_p = 0.05
    except:
        select_p = 0.05

    # display the top 5 destructive / protective phenoytpes
    df_variant_p = df_variant[df_variant["p-value"] <= select_p]
    df_variant_p_des = df_variant_p[df_variant_p["odds-ratio"] >= 1]
    df_variant_p_des = df_variant_p_des[["phewas phenotype", "gene_name", "snp", "odds-ratio", "p-value", "gwas-associations"]].reset_index().drop("index", axis= 1)
    df_variant_p_pro = df_variant_p[df_variant_p["odds-ratio"] < 1]
    df_variant_p_pro = df_variant_p_pro[["phewas phenotype", "gene_name", "snp", "odds-ratio", "p-value", "gwas-associations"]].reset_index().drop("index", axis= 1)

    st.header("Gene: " + "*" + select_gene + "*" + ", Variant: " + "*" + select_variant + "*" + ", P-value <= " + "*" + str(round(select_p, 4)) + "*")
    
    with st.container():
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Odds ratios of associated phenotypes")
            df_variant_p.sort_values(by=['odds-ratio'], inplace=True, ascending=True)
            fig = px.bar(df_variant_p, y = "phewas phenotype", x = "odds-ratio", color = "odds-ratio", color_continuous_scale = px.colors.sequential.RdBu_r, color_continuous_midpoint = 1, height= 1000)
            #fig.add_vline(x = 1)
            st.plotly_chart(fig, use_container_width= True)
        with col2:
            st.subheader("Data")
            st.markdown(get_table_download_link(df_variant_p), unsafe_allow_html=True)
            st.write('Risk allele-associated phenotypes (odds ratio > 1)')
            st.dataframe(df_variant_p_des, height = 400)
            st.write('Protective allele-associated phenotypes (odds ratio < 1)')
            st.dataframe(df_variant_p_pro, height = 400)

elif select == "Variant":
    st.markdown("This dashboard shows the associated phenotypes of your gene variants of interest.")

    # sidebar -- variant select box
    variant = sorted(df_selected["snp"].unique().tolist())
    select_variant = st.sidebar.selectbox('Variant', variant, key='4')

    # subset the data frame
    df_variant = df_selected[df_selected["snp"] == select_variant]

    # sidebar --  p-value slider
    df_variant.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
    select_p = st.sidebar.text_input(label = "P-value", help = "Defaults to p = 0.05. Accepts scientific notation, e.g., 5E-4, 3e-9", value = "0.05")
    try:
        if (float(select_p) <= 1) & (float(select_p) > 0):
            select_p = float(select_p)
        else:
            select_p = 0.05
    except:
        select_p = 0.05

    # display the top 5 destructive / protective phenoytpes
    df_variant_p = df_variant[df_variant["p-value"] <= select_p]
    df_variant_p_des = df_variant_p[df_variant_p["odds-ratio"] >= 1]
    df_variant_p_des = df_variant_p_des[["phewas phenotype", "gene_name", "snp", "odds-ratio", "p-value", "gwas-associations"]].reset_index().drop("index", axis= 1)
    df_variant_p_pro = df_variant_p[df_variant_p["odds-ratio"] < 1]
    df_variant_p_pro = df_variant_p_pro[["phewas phenotype", "gene_name", "snp", "odds-ratio", "p-value", "gwas-associations"]].reset_index().drop("index", axis= 1)

    st.header("Variant: " + "*" + select_variant + "*" + ", P-value <= " + "*" + str(round(select_p, 4)) + "*")
    st.write("Associated gene: " + ", ".join(df_variant_p["gene_name"].unique().tolist()))

    with st.container():
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Odds ratios of associated phenotypes")
            df_variant_p.sort_values(by=['odds-ratio'], inplace=True, ascending=True)
            fig = px.bar(df_variant_p, y = "phewas phenotype", x = "odds-ratio", color = "odds-ratio", color_continuous_scale = px.colors.sequential.RdBu_r, color_continuous_midpoint = 1, height= 1000)
            st.plotly_chart(fig, use_container_width= True)
        with col2:
            st.subheader("Data")
            st.markdown(get_table_download_link(df_variant_p), unsafe_allow_html=True)
            st.write('Risk allele-associated phenotypes (odds ratio > 1)')
            st.dataframe(df_variant_p_des, height = 400)
            st.write('Protective allele-associated phenotypes (odds ratio < 1)')
            st.dataframe(df_variant_p_pro, height = 400)

elif select == "GWAS":
    st.markdown("This dashboard shows the gene variants found in GWASs that are associated with your diseases of interest. It displays all associations, before separating associations into risk and protective. ")

    # extract diseases
    disease = ["--"]
    for i in list(set(df_selected["gwas-associations"].tolist())):
        disease.extend(i.split(", "))
    disease = sorted(list(set(disease)))

    # sidebar -- disease select box
    select_disease = st.sidebar.selectbox('Disease', disease, key='5')

    # subset the data frame for GWAS-associations
    df_disease_gwas = df_selected[df_selected["gwas-associations"].str.contains(select_disease)]
    df_disease_gwas = df_disease_gwas[["snp", "gene_name", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations"]]

    # sidebar --  p-value slider
    select_p = st.sidebar.text_input(label = "P-value", help = "Defaults to p = 0.05. Accepts scientific notation, e.g., 5E-4, 3e-9", value = "0.05")
    try:
        if (float(select_p) <= 1) & (float(select_p) > 0):
            select_p = float(select_p)
        else:
            select_p = 0.05
    except:
        select_p = 0.05
    df_disease_gwas = df_disease_gwas[df_disease_gwas["p-value"] <= select_p]
    
    # find the druggable genes
    # druggable evidence
    df_disease_gwas = pd.merge(df_disease_gwas, df_druggable, left_on='gene_name', right_on = "sym", how='left')
    df_disease_gwas['tdl'] = df_disease_gwas['tdl'].fillna("None")
    df_disease_gwas = df_disease_gwas.rename(columns={'tdl': 'druggability level'})

    # subset the data by odds ratio
    df_disease_gwas = df_disease_gwas.reset_index().drop("index", axis= 1)
    df_disease_gwas_sub_des = df_disease_gwas[df_disease_gwas["odds-ratio"] >= 1]
    df_disease_gwas_sub_pro = df_disease_gwas[df_disease_gwas["odds-ratio"] < 1]

    # count the druggability levels and generate data for pie chart (risk)
    df_tdl_des = df_disease_gwas_sub_des[["gene_name", "druggability level"]].drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
    df_tdl_pro = df_disease_gwas_sub_pro[["gene_name", "druggability level"]].drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
    df_tdl_all = df_disease_gwas[["gene_name", "druggability level"]].drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)

    st.header("Disease: " + "*" + select_disease + "*" + ", P-value <= " + "*" + str(round(select_p, 4)) + "*")

    st.subheader("All alleles")
    with st.container():
        col1, col2= st.columns([1,1])
        with col1:
            st.write("Percentage of genes in each druggability level:")
            elements_count = collections.Counter(df_tdl_all["druggability level"].tolist())
            gene_count = collections.Counter(df_tdl_all["gene_name"].tolist())
            df_disease_phewas_gwas_pie = pd.DataFrame(dict(elements_count).items(), columns= ["druggability level", "value"])
            gene_count_phewas_gwas_pie = pd.DataFrame(dict(gene_count).items(), columns= ["gene", "value"])
            labels = ["None", "Tdark", "Tbio", "Tchem", "Tclin"]
            df_disease_phewas_gwas_pie["druggability level"] = pd.Categorical(df_disease_phewas_gwas_pie["druggability level"], labels)
            fig = px.pie(df_disease_phewas_gwas_pie.sort_values("druggability level").reset_index(drop = True), values='value', names='druggability level', color='druggability level', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"})
            fig.update_traces(textposition='inside', textinfo='percent+label', showlegend=False, sort=False, rotation=0)
            st.plotly_chart(fig)
        with col2:
            fig = px.scatter(df_disease_gwas, x = "gene_name", y = "odds-ratio", color = "druggability level", hover_name='snp', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"},
                             category_orders = {
                    "druggability level": ["None", "Tdark", "Tbio", "Tchem", "Tclin"]})
            # size = [20]*len(df_disease_phewas_gwas),
            fig.add_hline(y = 1, line_width=1)
            st.plotly_chart(fig, use_container_width= True)
    with st.expander("See data"):
            st.write("GWAS data:")
            df_disease_gwas = df_disease_gwas[["gene_name", "snp", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations", "druggability level"]]
            df_disease_gwas.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
            df_disease_gwas = df_disease_gwas.reset_index().drop("index", axis= 1)
            st.markdown(get_table_download_link(df_disease_gwas.reset_index()), unsafe_allow_html=True)
            st.write(df_disease_gwas)

    st.subheader("Risk alleles (odds ratio > 1)")
    with st.container():
        col1, col2= st.columns([1,1])
        with col1:
            st.write("Percentage of genes in each druggability level:")
            elements_count = collections.Counter(df_tdl_des["druggability level"].tolist())
            df_disease_phewas_gwas_sub_des_pie = pd.DataFrame(dict(elements_count).items(), columns= ["druggability level", "value"])
            labels = ["None", "Tdark", "Tbio", "Tchem", "Tclin"]
            df_disease_phewas_gwas_sub_des_pie["druggability level"] = pd.Categorical(df_disease_phewas_gwas_sub_des_pie["druggability level"], labels)
            fig = px.pie(df_disease_phewas_gwas_sub_des_pie.sort_values("druggability level").reset_index(drop = True), values='value', names='druggability level', color='druggability level', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"})
            fig.update_traces(textposition='inside', textinfo='percent+label', showlegend=False, sort=False, rotation=0)
            st.plotly_chart(fig)
        with col2:
            fig = px.scatter(df_disease_gwas_sub_des, x = "gene_name", y = "odds-ratio", color = "druggability level", hover_name='snp', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"},
                             category_orders = {
                    "druggability level": ["None", "Tdark", "Tbio", "Tchem", "Tclin"]})
            # size = [20]*len(df_disease_phewas_gwas),
            fig.add_hline(y = 1, line_width=1)
            st.plotly_chart(fig, use_container_width= True)
    with st.expander("See data"):
            st.write("GWAS data (Odds-ratio > 1):")
            df_disease_gwas_sub_des = df_disease_gwas_sub_des[["gene_name", "snp", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations", "druggability level"]]
            df_disease_gwas_sub_des.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
            df_disease_gwas_sub_des = df_disease_gwas_sub_des.reset_index().drop("index", axis= 1)
            st.markdown(get_table_download_link(df_disease_gwas_sub_des.reset_index()), unsafe_allow_html=True)
            st.write(df_disease_gwas_sub_des)
    
    st.subheader("Protective alleles (odds ratio < 1)")
    with st.container():
        col1, col2= st.columns([1,1])
        with col1:
            st.write("Percentage of genes in each druggability level:")
            elements_count = collections.Counter(df_tdl_pro["druggability level"].tolist())
            df_disease_phewas_gwas_sub_pro_pie = pd.DataFrame(dict(elements_count).items(), columns= ["druggability level", "value"])
            labels = ["None", "Tdark", "Tbio", "Tchem", "Tclin"]
            df_disease_phewas_gwas_sub_pro_pie["druggability level"] = pd.Categorical(df_disease_phewas_gwas_sub_pro_pie["druggability level"], labels)
            fig = px.pie(df_disease_phewas_gwas_sub_pro_pie.sort_values("druggability level").reset_index(drop = True), values='value', names='druggability level', color='druggability level', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"})
            fig.update_traces(textposition='inside', textinfo='percent+label', showlegend=False, sort=False, rotation=0)
            st.plotly_chart(fig)
        with col2:
            fig = px.scatter(df_disease_gwas_sub_pro, x = "gene_name", y = "odds-ratio", color = "druggability level", hover_name='snp', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"},
                             category_orders = {
                    "druggability level": ["None", "Tdark", "Tbio", "Tchem", "Tclin"]})
            # size = [20]*len(df_disease_phewas_gwas),
            fig.add_hline(y = 1, line_width=1)
            st.plotly_chart(fig, use_container_width= True)
    with st.expander("See data"):
            st.write("GWAS data (Odds-ratio < 1):")
            df_disease_gwas_sub_pro = df_disease_gwas_sub_pro[["gene_name", "snp", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations", "druggability level"]]
            df_disease_gwas_sub_pro.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
            df_disease_gwas_sub_pro = df_disease_gwas_sub_pro.reset_index().drop("index", axis= 1)
            st.markdown(get_table_download_link(df_disease_gwas_sub_pro.reset_index()), unsafe_allow_html=True)
            st.write(df_disease_gwas_sub_pro)
            
elif select == "PheWAS":
    st.markdown("This dashboard shows the gene variants found in PheWASs that are associated with your diseases of interest. It displays all associations, before separating associations into risk and protective.")

    # extract diseases
    disease = ["--"]
    disease.extend(df_selected["phewas phenotype"].tolist())
    disease = sorted(list(set(disease)))

    # sidebar -- disease select box
    select_disease = st.sidebar.selectbox('Disease', disease, key='6')

    # subset the data frame for PheWAS
    df_disease_phewas = df_selected[df_selected["phewas phenotype"].str.contains(select_disease)]

    # subset the data frame for phewas phenotype
    df_disease_phewas = df_disease_phewas[["snp", "gene_name", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations"]]

    # sidebar --  p-value slider
    select_p = st.sidebar.text_input(label = "P-value", help = "Defaults to p = 0.05. Accepts scientific notation, e.g., 5E-4, 3e-9", value = "0.05")
    try:
        if (float(select_p) <= 1) & (float(select_p) > 0):
            select_p = float(select_p)
        else:
            select_p = 0.05
    except:
        select_p = 0.05
    df_disease_phewas = df_disease_phewas[df_disease_phewas["p-value"] <= select_p]
    
    # find the druggable genes
    # druggable evidence
    df_disease_phewas = pd.merge(df_disease_phewas, df_druggable, left_on='gene_name', right_on = "sym", how='left')
    df_disease_phewas['tdl'] = df_disease_phewas['tdl'].fillna("None")
    df_disease_phewas = df_disease_phewas.rename(columns={'tdl': 'druggability level'})

    # subset the data by odds ratio
    df_disease_phewas = df_disease_phewas.reset_index().drop("index", axis= 1)
    df_disease_phewas_sub_des = df_disease_phewas[df_disease_phewas["odds-ratio"] >= 1]
    df_disease_phewas_sub_pro = df_disease_phewas[df_disease_phewas["odds-ratio"] < 1]

    # count the druggability levels and generate data for pie chart (risk)
    df_tdl_des = df_disease_phewas_sub_des[["gene_name", "druggability level"]].drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
    df_tdl_pro = df_disease_phewas_sub_pro[["gene_name", "druggability level"]].drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
    df_tdl_all = df_disease_phewas[["gene_name", "druggability level"]].drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)

    st.header("Disease: " + "*" + select_disease + "*" + ", P-value <= " + "*" + str(round(select_p, 4)) + "*")

    st.subheader("All alleles")
    with st.container():
        col1, col2= st.columns([1,1])
        with col1:
            st.write("Percentage of genes in each druggability level:")
            elements_count = collections.Counter(df_tdl_all["druggability level"].tolist())
            df_disease_phewas_gwas_pie = pd.DataFrame(dict(elements_count).items(), columns= ["druggability level", "value"])
            labels = ["None", "Tdark", "Tbio", "Tchem", "Tclin"]
            df_disease_phewas_gwas_pie["druggability level"] = pd.Categorical(df_disease_phewas_gwas_pie["druggability level"], labels)
            fig = px.pie(df_disease_phewas_gwas_pie.sort_values("druggability level").reset_index(drop = True), values='value', names='druggability level', color='druggability level', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"})
            fig.update_traces(textposition='inside', textinfo='percent+label', showlegend=False, sort=False, rotation=0)
            st.plotly_chart(fig)
        with col2:
            fig = px.scatter(df_disease_phewas, x = "gene_name", y = "odds-ratio", color = "druggability level", hover_name='snp', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"},
                             category_orders = {
                    "druggability level": ["None", "Tdark", "Tbio", "Tchem", "Tclin"]})
            # size = [20]*len(df_disease_phewas_gwas),
            fig.add_hline(y = 1, line_width=1)
            st.plotly_chart(fig, use_container_width= True)
    with st.expander("See data"):
            st.write("PheWAS data:")
            df_disease_phewas = df_disease_phewas[["gene_name", "snp", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations", "druggability level"]]
            df_disease_phewas.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
            df_disease_phewas = df_disease_phewas.reset_index().drop("index", axis= 1)
            st.markdown(get_table_download_link(df_disease_phewas.reset_index()), unsafe_allow_html=True)
            st.write(df_disease_phewas)

    st.subheader("Risk alleles (odds ratio > 1)")
    with st.container():
        col1, col2= st.columns([1,1])
        with col1:
            st.write("Percentage of genes in each druggability level:")
            elements_count = collections.Counter(df_tdl_des["druggability level"].tolist())
            df_disease_phewas_gwas_sub_des_pie = pd.DataFrame(dict(elements_count).items(), columns= ["druggability level", "value"])
            labels = ["None", "Tdark", "Tbio", "Tchem", "Tclin"]
            df_disease_phewas_gwas_sub_des_pie["druggability level"] = pd.Categorical(df_disease_phewas_gwas_sub_des_pie["druggability level"], labels)
            fig = px.pie(df_disease_phewas_gwas_sub_des_pie.sort_values("druggability level").reset_index(drop = True), values='value', names='druggability level', color='druggability level', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"})
            fig.update_traces(textposition='inside', textinfo='percent+label', showlegend=False, sort=False, rotation=0)
            st.plotly_chart(fig)
        with col2:
            fig = px.scatter(df_disease_phewas_sub_des, x = "gene_name", y = "odds-ratio", color = "druggability level", hover_name='snp', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"},
                             category_orders = {
                    "druggability level": ["None", "Tdark", "Tbio", "Tchem", "Tclin"]})
            # size = [20]*len(df_disease_phewas_gwas),
            fig.add_hline(y = 1, line_width=1)
            st.plotly_chart(fig, use_container_width= True)
    with st.expander("See data"):
            st.write("PheWAS data (Odds-ratio >= 1):")
            df_disease_phewas_sub_des = df_disease_phewas_sub_des[["gene_name", "snp", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations", "druggability level"]]
            df_disease_phewas_sub_des.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
            df_disease_phewas_sub_des = df_disease_phewas_sub_des.reset_index().drop("index", axis= 1)
            st.markdown(get_table_download_link(df_disease_phewas_sub_des.reset_index()), unsafe_allow_html=True)
            st.write(df_disease_phewas_sub_des)
    
    st.subheader("Protective alleles (odds ratio < 1)")
    with st.container():
        col1, col2= st.columns([1,1])
        with col1:
            st.write("Percentage of genes in each druggability level:")
            elements_count = collections.Counter(df_tdl_pro["druggability level"].tolist())
            df_disease_phewas_gwas_sub_pro_pie = pd.DataFrame(dict(elements_count).items(), columns= ["druggability level", "value"])
            labels = ["None", "Tdark", "Tbio", "Tchem", "Tclin"]
            df_disease_phewas_gwas_sub_pro_pie["druggability level"] = pd.Categorical(df_disease_phewas_gwas_sub_pro_pie["druggability level"], labels)
            fig = px.pie(df_disease_phewas_gwas_sub_pro_pie.sort_values("druggability level").reset_index(drop = True), values='value', names='druggability level', color='druggability level', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"})
            fig.update_traces(textposition='inside', textinfo='percent+label', showlegend=False, sort=False, rotation=0)
            st.plotly_chart(fig)
        with col2:
            fig = px.scatter(df_disease_phewas_sub_pro, x = "gene_name", y = "odds-ratio", color = "druggability level", hover_name='snp', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"},
                             category_orders = {
                    "druggability level": ["None", "Tdark", "Tbio", "Tchem", "Tclin"]})
            # size = [20]*len(df_disease_phewas_gwas),
            fig.add_hline(y = 1, line_width=1)
            st.plotly_chart(fig, use_container_width= True)
    with st.expander("See data"):
            st.write("PheWAS data (Odds-ratio < 1):")
            df_disease_phewas_sub_pro = df_disease_phewas_sub_pro[["gene_name", "snp", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations", "druggability level"]]
            df_disease_phewas_sub_pro.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
            df_disease_phewas_sub_pro = df_disease_phewas_sub_pro.reset_index().drop("index", axis= 1)
            st.markdown(get_table_download_link(df_disease_phewas_sub_pro.reset_index()), unsafe_allow_html=True)
            st.write(df_disease_phewas_sub_pro)

elif select == "GWAS_PheWAS Union":
    st.markdown("This dashboard shows the gene variants found in either GWASs or PheWASs or both that are associated with your diseases of interest. It displays all associations, before separating associations into risk and protective.")
    
    # extract diseases
    disease = ["--"]
    for i in list(set(df_selected["gwas-associations"].tolist())):
        disease.extend(i.split(", "))
    disease.extend(df_selected["phewas phenotype"].tolist())
    disease = sorted(list(set(disease)))

    # sidebar -- disease select box
    select_disease = st.sidebar.selectbox('Disease', disease, key='7')

    # subset the dataframe with the aggregation of phewas and gwas association
    df_disease_phewas_gwas = df_selected[df_selected["phewas phenotype"].str.contains(select_disease) | df_selected["gwas-associations"].str.contains(select_disease)]

    # subset the data frame for phewas phenotype and gwas association
    df_disease_phewas_gwas = df_disease_phewas_gwas[["gene_name", "snp", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations"]]

    # sidebar --  p-value slider
    select_p = st.sidebar.text_input(label = "P-value", help = "Defaults to p = 0.05. Accepts scientific notation, e.g., 5E-4, 3e-9", value = "0.05")
    try:
        if (float(select_p) <= 1) & (float(select_p) > 0):
            select_p = float(select_p)
        else:
            select_p = 0.05
    except:
        select_p = 0.05
    df_disease_phewas_gwas = df_disease_phewas_gwas[df_disease_phewas_gwas["p-value"] <= select_p]
    
    # look for evidence of druggability    
    # druggable evidence
    df_disease_phewas_gwas = pd.merge(df_disease_phewas_gwas, df_druggable, left_on='gene_name', right_on = "sym", how='left')
    df_disease_phewas_gwas['tdl'] = df_disease_phewas_gwas['tdl'].fillna("None")
    df_disease_phewas_gwas = df_disease_phewas_gwas.rename(columns={'tdl': 'druggability level'})

    # subset the data by odds ratio
    df_disease_phewas_gwas = df_disease_phewas_gwas.reset_index().drop("index", axis= 1)
    df_disease_phewas_gwas_sub_des = df_disease_phewas_gwas[df_disease_phewas_gwas["odds-ratio"] >= 1]
    df_disease_phewas_gwas_sub_pro = df_disease_phewas_gwas[df_disease_phewas_gwas["odds-ratio"] < 1]

    # count the druggability levels and generate data for pie chart (risk)
    df_tdl_des = df_disease_phewas_gwas_sub_des[["gene_name", "druggability level"]].drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
    df_tdl_pro = df_disease_phewas_gwas_sub_pro[["gene_name", "druggability level"]].drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
    df_tdl_all = df_disease_phewas_gwas[["gene_name", "druggability level"]].drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)

    st.header("Disease: " + "*" + select_disease + "*" + ", P-value <= " + "*" + str(round(select_p, 4)) + "*")

    st.subheader("All alleles")
    with st.container():
        col1, col2= st.columns([1,1])
        with col1:
            st.write("Percentage of genes in each druggability level:")
            elements_count = collections.Counter(df_tdl_all["druggability level"].tolist())
            df_disease_phewas_gwas_pie = pd.DataFrame(dict(elements_count).items(), columns= ["druggability level", "value"])
            labels = ["None", "Tdark", "Tbio", "Tchem", "Tclin"]
            df_disease_phewas_gwas_pie["druggability level"] = pd.Categorical(df_disease_phewas_gwas_pie["druggability level"], labels)
            fig = px.pie(df_disease_phewas_gwas_pie.sort_values("druggability level").reset_index(drop = True), values='value', names='druggability level', color='druggability level', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"})
            fig.update_traces(textposition='inside', textinfo='percent+label', showlegend=False, sort=False, rotation=0)
            st.plotly_chart(fig)
        with col2:
            fig = px.scatter(df_disease_phewas_gwas, x = "gene_name", y = "odds-ratio", color = "druggability level", hover_name='snp', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"},
                             category_orders = {
                    "druggability level": ["None", "Tdark", "Tbio", "Tchem", "Tclin"]})
            # size = [20]*len(df_disease_phewas_gwas),
            fig.add_hline(y = 1, line_width=1)
            st.plotly_chart(fig, use_container_width= True)
    with st.expander("See data"):
            st.write("PheWAS+GWAS data:")
            df_disease_phewas_gwas = df_disease_phewas_gwas[["gene_name", "snp", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations", "druggability level"]]
            df_disease_phewas_gwas.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
            df_disease_phewas_gwas = df_disease_phewas_gwas.reset_index().drop("index", axis= 1)
            st.markdown(get_table_download_link(df_disease_phewas_gwas.reset_index()), unsafe_allow_html=True)
            st.write(df_disease_phewas_gwas)

    st.subheader("Risk alleles (odds ratio > 1)")
    with st.container():
        col1, col2= st.columns([1,1])
        with col1:
            st.write("Percentage of genes in each druggability level:")
            elements_count = collections.Counter(df_tdl_des["druggability level"].tolist())
            df_disease_phewas_gwas_sub_des_pie = pd.DataFrame(dict(elements_count).items(), columns= ["druggability level", "value"])
            labels = ["None", "Tdark", "Tbio", "Tchem", "Tclin"]
            df_disease_phewas_gwas_sub_des_pie["druggability level"] = pd.Categorical(df_disease_phewas_gwas_sub_des_pie["druggability level"], labels)
            fig = px.pie(df_disease_phewas_gwas_sub_des_pie.sort_values("druggability level").reset_index(drop = True), values='value', names='druggability level', color='druggability level', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"})
            fig.update_traces(textposition='inside', textinfo='percent+label', showlegend=False, sort=False, rotation=0)
            st.plotly_chart(fig)
        with col2:
            fig = px.scatter(df_disease_phewas_gwas_sub_des, x = "gene_name", y = "odds-ratio", color = "druggability level", hover_name='snp', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"},
                             category_orders = {
                    "druggability level": ["None", "Tdark", "Tbio", "Tchem", "Tclin"]})
            # size = [20]*len(df_disease_phewas_gwas),
            fig.add_hline(y = 1, line_width=1)
            st.plotly_chart(fig, use_container_width= True)
    with st.expander("See data"):
            st.write("PheWAS+GWAS data (Odds-ratio >= 1):")
            df_disease_phewas_gwas_sub_des = df_disease_phewas_gwas_sub_des[["gene_name", "snp", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations", "druggability level"]]
            df_disease_phewas_gwas_sub_des.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
            df_disease_phewas_gwas_sub_des = df_disease_phewas_gwas_sub_des.reset_index().drop("index", axis= 1)
            st.markdown(get_table_download_link(df_disease_phewas_gwas_sub_des.reset_index()), unsafe_allow_html=True)
            st.write(df_disease_phewas_gwas_sub_des)
    
    st.subheader("Protective alleles (odds ratio < 1)")
    with st.container():
        col1, col2= st.columns([1,1])
        with col1:
            st.write("Percentage of genes in each druggability level:")
            elements_count = collections.Counter(df_tdl_pro["druggability level"].tolist())
            df_disease_phewas_gwas_sub_pro_pie = pd.DataFrame(dict(elements_count).items(), columns= ["druggability level", "value"])
            labels = ["None", "Tdark", "Tbio", "Tchem", "Tclin"]
            df_disease_phewas_gwas_sub_pro_pie["druggability level"] = pd.Categorical(df_disease_phewas_gwas_sub_pro_pie["druggability level"], labels)
            fig = px.pie(df_disease_phewas_gwas_sub_pro_pie.sort_values("druggability level").reset_index(drop = True), values='value', names='druggability level', color='druggability level', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"})
            fig.update_traces(textposition='inside', textinfo='percent+label', showlegend=False, sort=False, rotation=0)
            st.plotly_chart(fig)
        with col2:
            fig = px.scatter(df_disease_phewas_gwas_sub_pro, x = "gene_name", y = "odds-ratio", color = "druggability level", hover_name='snp', color_discrete_map = {
                "None": "rgb(203,213,232)",
                "Tdark": "rgb(141,160,203)",
                "Tbio": "rgb(223,217,164)",
                "Tchem": "rgb(229,134,6)",
                "Tclin": "#DC3912"},
                             category_orders = {
                    "druggability level": ["None", "Tdark", "Tbio", "Tchem", "Tclin"]})
            # size = [20]*len(df_disease_phewas_gwas),
            fig.add_hline(y = 1, line_width=1)
            st.plotly_chart(fig, use_container_width= True)
    with st.expander("See data"):
            st.write("PheWAS+GWAS data (Odds-ratio < 1):")
            df_disease_phewas_gwas_sub_pro = df_disease_phewas_gwas_sub_pro[["gene_name", "snp", "odds-ratio", "p-value", "phewas phenotype", "gwas-associations", "druggability level"]]
            df_disease_phewas_gwas_sub_pro.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
            df_disease_phewas_gwas_sub_pro = df_disease_phewas_gwas_sub_pro.reset_index().drop("index", axis= 1)
            st.markdown(get_table_download_link(df_disease_phewas_gwas_sub_pro.reset_index()), unsafe_allow_html=True)
            st.write(df_disease_phewas_gwas_sub_pro)
    
elif select == "GWAS_PheWAS Intersection":
    st.markdown("This dashboard shows the gene variants found in **both** GWASs and PheWASs that are associated with your diseases of interest. Gene variants that lie in this intersection are then further analysed in their druggability, association odds-ratio, protein-protein interactions and gene ontology term enrichment.")

    # extract diseases
    disease = ["--"]
    for i in list(set(df_selected["gwas-associations"].tolist())):
        disease.extend(i.split(", "))
    phewas_disease_list = list(set(df_selected["phewas phenotype"].tolist()))
    intersection_list = sorted(list(set(disease) & set(phewas_disease_list)))
    
    # sidebar -- disease select box
    select_disease = st.sidebar.selectbox('Disease', intersection_list, key='8')

    # sidebar --  p-value slider
    select_p = st.sidebar.text_input(label = "P-value", help = "Defaults to p = 0.05. Accepts scientific notation, e.g., 5E-4, 3e-9", value = "0.05")
    try:
        if (float(select_p) <= 1) & (float(select_p) > 0):
            select_p = float(select_p)
        else:
            select_p = 0.05
    except:
        select_p = 0.05
    df_selected = df_selected[df_selected["p-value"] <= select_p]

    # subset the data frame for either PheWAS and GWAS
    df_disease_phewas_or_gwas = df_selected[df_selected["phewas phenotype"].str.contains(select_disease) | df_selected["gwas-associations"].str.contains(select_disease)]
    df_disease_phewas_or_gwas.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
    df_disease_phewas_or_gwas = df_disease_phewas_or_gwas[["snp", "gene_name", "phewas phenotype", "odds-ratio", "p-value", "gwas-associations"]].reset_index().drop("index", axis= 1)

    # subset the data frame for PheWAS
    df_disease_phewas = df_selected[df_selected["phewas phenotype"].str.contains(select_disease)]
    df_disease_phewas = df_disease_phewas.reset_index().drop("index", axis= 1)
    #df_disease_phewas["gene_snp"] = df_disease_phewas["gene_name"] + " " + df_disease_phewas['snp']
    df_disease_phewas.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
    df_disease_phewas_sub = df_disease_phewas[["snp", "gene_name", "phewas phenotype", "odds-ratio", "p-value", "gwas-associations"]].reset_index().drop("index", axis= 1)

    # subset the data frame for GWAS-associations
    df_disease_gwas = df_selected[df_selected["gwas-associations"].str.contains(select_disease)]
    df_disease_gwas = df_disease_gwas.reset_index().drop("index", axis= 1)
    #df_disease_gwas["gene_snp"] = df_disease_gwas["gene_name"] + " " + df_disease_gwas['snp']
    df_disease_gwas.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
    df_disease_gwas_sub = df_disease_gwas[["snp", "gene_name", "phewas phenotype", "odds-ratio", "p-value", "gwas-associations"]].reset_index().drop("index", axis= 1)

    # subset the data frame with overlapped phewas and gwas association
    df_disease_phewas_gwas = df_selected[df_selected["phewas phenotype"].str.contains(select_disease) & df_selected["gwas-associations"].str.contains(select_disease)]
    df_disease_phewas_gwas = df_disease_phewas_gwas.reset_index().drop("index", axis= 1)
    #df_disease_phewas_gwas["gene_snp"] = df_disease_phewas_gwas["gene_name"] + " " + df_disease_phewas_gwas['snp']
    df_disease_phewas_gwas.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
    df_disease_phewas_gwas_sub = df_disease_phewas_gwas[["snp", "gene_name", "phewas phenotype", "odds-ratio", "p-value", "gwas-associations"]].reset_index().drop("index", axis= 1)

    # find the druggable genes
    druggable_gene = []
    for g in df_disease_phewas_gwas["gene_name"].tolist():
        if g in df_druggable["sym"].tolist():
            druggable_gene.append(g)
    try:
        df_disease_phewas_gwas_druggable = df_disease_phewas_gwas[df_disease_phewas_gwas['gene_name'].isin(druggable_gene)].reset_index().drop("index", axis= 1)
        #df_disease_phewas_gwas_druggable["gene_snp"] = df_disease_phewas_gwas_druggable["gene_name"] + " " + df_disease_phewas_gwas_druggable['snp']
        df_disease_phewas_gwas_druggable = pd.merge(df_disease_phewas_gwas_druggable,df_druggable, left_on='gene_name', right_on = "sym", how='left')
        df_disease_phewas_gwas_druggable.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
        df_disease_phewas_gwas_druggable_sub = df_disease_phewas_gwas_druggable[["gene_name", "snp", "phewas phenotype", "odds-ratio", "p-value", "gwas-associations", "tdl"]].reset_index().drop("index", axis= 1)
        df_disease_phewas_gwas_druggable_sub = df_disease_phewas_gwas_druggable_sub.rename(columns={'tdl': 'druggability level'})
        druggable_rate = len(df_disease_phewas_gwas_druggable_sub)/len(df_disease_phewas_gwas_sub)
        non_druggable_count = len(df_disease_phewas_gwas_sub) - len(df_disease_phewas_gwas_druggable_sub)
        non_druggable_rate = non_druggable_count/len(df_disease_phewas_gwas_sub)

        st.header("Disease: " + "*" + select_disease + "*" + ", P-value <= " + "*" + str(round(select_p, 4)) + "*")
        
        # get the phewas and gwas data
        phewas_set = df_disease_phewas["snp"].unique().tolist()
        gwas_set = df_disease_gwas["snp"].unique().tolist()
        phewas_only = np.setdiff1d(phewas_set,gwas_set)
        gwas_only = np.setdiff1d(gwas_set,phewas_set)
        phewas_gwas = set(phewas_set) & set(gwas_set)

        option_phewas_gwas_pie = {
            "tooltip": {
                "trigger": 'item',
                "formatter": '{b} : {c} ({d}%)',
            },
            "legend": {
                "orient": 'vertical',
                "left": 'left'
            },
            "grid": { "containLabel": True },
            "series": [
                {
                    "type": 'pie',
                    "radius": 150,
                    "center": ['50%', '50%'],
                    "selectedMode": 'single',
                    "clockwise": True,
                    "labelLayout":{
                        "draggable": True
                    },
                    "labelLine": {
                        "length": 30,
                    },
                    "itemStyle": {
                        "borderWidth": 3
                    },
                    "data": [
                        {
                            "value": len(phewas_gwas),
                            "name": 'PheWAS_GWAS intersection',
                            "selected": True,
                            "label": {
                                "formatter": '\n'.join([
                                    '{title|{b}}{abg|}',
                                    '  {druggabilityHead|Druggability}{valueHead|Count}{rateHead|Rate}',
                                    '{hr|}',
                                    '  {druggability|Druggable}{value|' + str(len(df_disease_phewas_gwas_druggable_sub)) + '}{rate|' + str(round(druggable_rate, 2)) +'}',
                                    '  {Nondruggability|Non-druggable}{value|' + str(non_druggable_count) + '}{rate|' + str(round(non_druggable_rate, 2)) +'}'
                                ]),
                                "backgroundColor": 'white',
                                "borderColor": '#777',
                                "borderWidth": 1,
                                "borderRadius": 4,
                                "rich": {
                                    "title": {
                                        "color": '#eee',
                                        "align": 'center'
                                    },
                                    "abg": {
                                        "backgroundColor": '#222',
                                        "width": '100%',
                                        "align": 'right',
                                        "height": 25,
                                        "borderRadius": [4, 4, 0, 0]
                                    },
                                    "druggability": {
                                        "height": 20,
                                        "width": 70,
                                        "align": 'center',
                                        "lineHeight": 30,
                                        "color": '#fff',
                                        "borderRadius": 5,
                                        "backgroundColor": '#4C5058',
                                        "padding": [0, 4],
                                    },
                                    "Nondruggability": {
                                        "height": 20,
                                        "width": 70,
                                        "align": 'center',
                                        "lineHeight": 20,
                                    },
                                    "druggabilityHead": {
                                        "color": '#333',
                                        "height": 24,
                                        "width": 10,
                                        "align": 'left'
                                    },
                                    "hr": {
                                        "borderColor": '#777',
                                        "width": '100%',
                                        "borderWidth": 0.5,
                                        "height": 0
                                    },
                                    "value": {
                                        "width": 10,
                                        "padding": [0, 20, 0, 30],
                                        "align": 'right'
                                    },
                                    "valueHead": {
                                        "color": '#333',
                                        "width": 10,
                                        "padding": [0, 20, 0, 30],
                                        "align": 'right'
                                    },
                                    "rate": {
                                        "width": 20,
                                        "align": 'right',
                                        "padding": [0, 10, 0, 0]
                                    },
                                    "rateHead": {
                                        "color": '#333',
                                        "width": 20,
                                        "align": 'right',
                                        "padding": [0, 10, 0, 0]
                                    }
                                },
                                "offset": [0, -10],
                            },
                            "labelLine":{
                                "showAbove":True,
                                "length": 100,
                                "length2": 50
                            }
                        },
                        {"value": len(phewas_only), "name": 'PheWAS only', "label":{
                            "show": True,
                            "formatter": '{b} : {c} ({d}%)'
                            }},
                        {"value": len(gwas_only), "name": 'GWAS only', "label":{
                            "show": True,
                            "formatter": '{b} : {c} ({d}%)'
                            }}
                    ],
                    "emphasis": {
                        "itemStyle": {
                            "shadowBlur": 10,
                            "shadowOffsetX": 0,
                            "shadowColor": 'rgba(0, 0, 0, 0.5)'
                        }
                    }
                , "color": ["#DC3912", "rgb(203,213,232)", "rgb(141,160,203)"]}
            ]
        }   
        
        # generate snp, gene, and druggable level data for pie chart
        df_druglvl = df_druggable[["sym", "tdl"]]
        df_druglvl = df_druglvl.drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
        
        # druggable evidence
        df_disease_phewas_gwas = pd.merge(df_disease_phewas_gwas, df_druglvl, left_on='gene_name', right_on = "sym", how='left')
        df_disease_phewas_gwas['tdl'] = df_disease_phewas_gwas['tdl'].fillna("None")
        df_disease_phewas_gwas = df_disease_phewas_gwas.rename(columns={'tdl': 'druggability level'})

        # get the druggable snp and gene data
        drug_level = df_disease_phewas_gwas["druggability level"].unique().tolist()
        druggable_data = []
        for d in drug_level:

            if d == "None":
                colour = "rgb(203,213,232)"
            elif d == "Tdark":
                colour = "rgb(141,160,203)"
            elif d == "Tbio":
                colour = "rgb(223,217,164)"
            elif d == "Tchem":
                colour = "rgb(229,134,6)"
            elif d == "Tclin":
                colour = "#DC3912"
            
            gene = df_disease_phewas_gwas[df_disease_phewas_gwas["druggability level"] == d]["gene_name"].unique().tolist()
            parent = []
            for i in gene:
                snp = df_disease_phewas_gwas[df_disease_phewas_gwas["gene_name"] == i]["snp"].unique().tolist()
                children = []
                for y in snp:
                    chil_info = {
                        "name": y,
                        "value": df_disease_phewas_gwas[df_disease_phewas_gwas["snp"] == y]["odds-ratio"].tolist()[0],
                        "itemStyle": {
                            "color": colour
                        }
                    }
                    children.append(chil_info)
                child_info = {
                    "name": i,
                    "children": children,
                    "itemStyle": {
                    "color": colour
                    }
                }
                parent.append(child_info)
            parant_infor = {
                "name": d,
                "children": parent,
                "itemStyle": {
                    "color": colour
                    }
            }
            druggable_data.append(parant_infor)
        
        option_druggable = {
            "title": {
                #"text": 'Gene Level and SNP Level',
                "textStyle": {
                    "fontSize": 14,
                    "align": 'center'
                },
                "subtextStyle": {
                    "align": 'center'
                },
            },
            "series": {
                "type": 'sunburst',

                "data": druggable_data,
                "radius": [20, 200],
                "center": ['50%', '50%'],
                "sort": "desc",

                "emphasis": {
                    "focus": 'descendant'
                },
                "animation": True,
                "levels": [{}, {
                    "r0": '15%',
                    "r": '35%',
                    "label": {
                        "rotate": 'radial'
                    },
                    "itemStyle": {
                        "borderWidth": 3
                    }
                }, {
                    "r0": '35%',
                    "r": '65%',
                    "label": {
                        "align": 'right'
                    },
                    "itemStyle": {
                        "borderWidth": 3
                    }
                }, {
                    "r0": '65%',
                    "r": '70%',
                    "label": {
                        "position": 'outside',
                        "padding": 0,
                        "silent": False
                    },
                    "itemStyle": {
                        "borderWidth": 3
                    }
                }]
            }
        }
        
        with st.container():
            col1, col2= st.columns(2)
            with col1:
                st.write("Percentage of SNPs from PheWASs and GWASs")
                st_echarts(options = option_phewas_gwas_pie, key= '2', height = 400)

            with col2:
                # druggable genes bar chart
                if not df_disease_phewas_gwas_druggable.empty:
                    # generate bar chart for druggable genes
                    st.write("Druggability level of SNPs from PheWAS & GWAS intersection")
                    st_echarts(options= option_druggable, key= '3', height = 400)

        snp_druggable = []
        for i in range(len(df_disease_phewas_gwas_druggable_sub["snp"])):
            snp_info = {
                "value": df_disease_phewas_gwas_druggable_sub["odds-ratio"].tolist()[i],
                "name": df_disease_phewas_gwas_druggable_sub["snp"].tolist()[i]
            }
            snp_druggable.append(snp_info)

        with st.container():
            col1, col2= st.columns([1,1])
            with col1:
                st.write("PheWAS & GWAS intersection alleles odds ratio:")
                fig = px.scatter(df_disease_phewas_gwas, x = "gene_name", y = "odds-ratio", color = "snp")
                fig.add_hline(y = 1, line_width=1)
                st.plotly_chart(fig, use_container_width= True)
            with col2:
                st.write("Protein-protein interaction network of PheWAS & GWAS intersection alleles:")
                g=net.Network(height='400px', width='700px')
                select_gene = df_disease_phewas_gwas["gene_name"].tolist()
                try:
                    df_protein_interaction = proteins_interaction(input_protein= select_gene)
                    for i in df_protein_interaction.index:
                        g.add_node(df_protein_interaction["preferredName_A"][i])
                        g.add_node(df_protein_interaction["preferredName_B"][i])
                        g.add_edge(df_protein_interaction["preferredName_A"][i],df_protein_interaction["preferredName_B"][i])
                except Exception:
                    pass
                pv_static(g)

        with st.expander("GO enrichment analysis"): 
            df_disease_phewas_gwas =df_disease_phewas_gwas[["gene_name", "snp", "phewas phenotype", "odds-ratio", "p-value", "gwas-associations", "druggability level"]]
            df_enrichment = go_enrichment(input_gene = df_disease_phewas_gwas_druggable_sub['gene_name'].tolist())
            st.write("GO enrichment analysis")
            st.markdown(get_table_download_link(df_enrichment.reset_index()), unsafe_allow_html=True)
            st.dataframe(df_enrichment)
            st.write("PheWAS GWAS intersect data ")
            st.markdown(get_table_download_link(df_disease_phewas_gwas.reset_index()), unsafe_allow_html=True)
            st.dataframe(df_disease_phewas_gwas)
    except:
        st.subheader("No data found. Please try selecting another disease!")


elif select == "Protein-protein Interaction":
    st.markdown("This dashboard shows the protein-protein interaction networks and gene ontology enrichment for your diseases of interest.")

    # extract diseases
    disease = ["--"]
    for i in list(set(df_selected["gwas-associations"].tolist())):
        disease.extend(i.split(", "))
    disease.extend(df_selected["phewas phenotype"].tolist())
    disease = sorted(list(set(disease)))

    # sidebar -- disease select box
    select_disease = st.sidebar.selectbox('Disease', disease, key='9')

    # sidebar --  p-value slider
    select_p = st.sidebar.text_input(label = "P-value", help = "Defaults to p = 0.05. Accepts scientific notation, e.g., 5E-4, 3e-9", value = "0.05")
    try:
        if (float(select_p) <= 1) & (float(select_p) > 0):
            select_p = float(select_p)
        else:
            select_p = 0.05
    except:
        select_p = 0.05
    df_selected = df_selected[df_selected["p-value"] <= select_p]

    # sidebar -- data
    select_study = st.sidebar.selectbox('Study', ["PheWAS + GWAS", "PheWAS", "GWAS"], key='10')

    if select_study == "PheWAS + GWAS":
    # subset the data frame for the disease appearing in either PheWAS and GWAS
        df_disease_phewas_or_gwas = df_selected[df_selected["phewas phenotype"].str.contains(select_disease) | df_selected["gwas-associations"].str.contains(select_disease)]
        df_disease_phewas_or_gwas.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
        df_disease_phewas_or_gwas = df_disease_phewas_or_gwas[["gene_name", "snp", "phewas phenotype", "odds-ratio", "p-value", "gwas-associations"]].reset_index().drop("index", axis= 1)
        df_disease_phewas_or_gwas = df_disease_phewas_or_gwas
    elif select_study == "PheWAS":
        df_disease_phewas_or_gwas = df_selected[df_selected["phewas phenotype"].str.contains(select_disease)]
    elif select_study == "GWAS":
        df_disease_phewas_or_gwas = df_selected[df_selected["gwas-associations"].str.contains(select_disease)]

    # subset the data by odds ratio
    df_disease_phewas_or_gwas = df_disease_phewas_or_gwas.reset_index().drop("index", axis= 1)
    df_disease_phewas_or_gwas_des = df_disease_phewas_or_gwas[df_disease_phewas_or_gwas["odds-ratio"] >= 1]
    df_disease_phewas_or_gwas_pro = df_disease_phewas_or_gwas[df_disease_phewas_or_gwas["odds-ratio"] < 1]

    
    # protein-protein network (des genes)
    gene_set = df_disease_phewas_or_gwas_des["gene_name"].unique().tolist()
    g_des=net.Network(width='650px')
    try:
        df_protein_interaction = proteins_interaction(input_protein= gene_set)
        
        # compute the degree of nodes
        frequency_a = collections.Counter(df_protein_interaction["preferredName_A"].tolist())
        dict_frequency_a = dict(frequency_a)
        frequency_b = collections.Counter(df_protein_interaction["preferredName_B"].tolist())
        dict_frequency_b = dict(frequency_b)

        gene_degree = {}
        for i in gene_set:
            degree_a = 0 if dict_frequency_a.get(i) is None else dict_frequency_a.get(i)
            degree_b = 0 if dict_frequency_b.get(i) is None else dict_frequency_b.get(i)
            gene_degree[i] = degree_a + degree_b
        df_gene_degree_des = pd.DataFrame(gene_degree.items(), columns=['gene','degree'])
        df_gene_degree_des.sort_values(by=['degree'], inplace=True, ascending=False)
        df_gene_degree_des = df_gene_degree_des.set_index('gene')

        # add degree score
        df_disease_phewas_or_gwas_des = pd.merge(df_disease_phewas_or_gwas_des, df_gene_degree_des, left_on='gene_name', right_on = "gene", how='left')

        # generate PPI network
        for i in df_protein_interaction.index:
            g_des.add_node(df_protein_interaction["preferredName_A"][i], color = "rgb(229,134,6)")
            g_des.add_node(df_protein_interaction["preferredName_B"][i], color = "rgb(229,134,6)")
            g_des.add_edge(df_protein_interaction["preferredName_A"][i],df_protein_interaction["preferredName_B"][i], weight= df_protein_interaction["score"][i])
    except Exception:
        pass

    # protein-protein network (pro genes)
    gene_set = df_disease_phewas_or_gwas_pro["gene_name"].unique().tolist()
    g_pro=net.Network(width='650px')
    try:
        df_protein_interaction = proteins_interaction(input_protein= gene_set)
        
        # compute the degree of nodes
        frequency_a = collections.Counter(df_protein_interaction["preferredName_A"].tolist())
        dict_frequency_a = dict(frequency_a)
        frequency_b = collections.Counter(df_protein_interaction["preferredName_B"].tolist())
        dict_frequency_b = dict(frequency_b)

        gene_degree = {}
        for i in gene_set:
            degree_a = 0 if dict_frequency_a.get(i) is None else dict_frequency_a.get(i)
            degree_b = 0 if dict_frequency_b.get(i) is None else dict_frequency_b.get(i)
            gene_degree[i] = degree_a + degree_b
        df_gene_degree_pro = pd.DataFrame(gene_degree.items(), columns=['gene','degree'])
        df_gene_degree_pro.sort_values(by=['degree'], inplace=True, ascending=False)
        df_gene_degree_pro = df_gene_degree_pro.set_index('gene')
        
        # add degree score
        df_disease_phewas_or_gwas_pro = pd.merge(df_disease_phewas_or_gwas_pro, df_gene_degree_pro, left_on='gene_name', right_on = "gene", how='left')

        # generate PPI network
        for i in df_protein_interaction.index:
            g_pro.add_node(df_protein_interaction["preferredName_A"][i], color = "rgb(148, 203, 141)")
            g_pro.add_node(df_protein_interaction["preferredName_B"][i], color = "rgb(148, 203, 141)")
            g_pro.add_edge(df_protein_interaction["preferredName_A"][i],df_protein_interaction["preferredName_B"][i], weight= df_protein_interaction["score"][i])
    except Exception:
        pass

    # protein-protein network (overall genes)
    gene_set = df_disease_phewas_or_gwas["gene_name"].unique().tolist()
    g_all=net.Network(width='650px')
    try:
        df_protein_interaction = proteins_interaction(input_protein= gene_set)
        
        # compute the degree of nodes
        frequency_a = collections.Counter(df_protein_interaction["preferredName_A"].tolist())
        dict_frequency_a = dict(frequency_a)
        frequency_b = collections.Counter(df_protein_interaction["preferredName_B"].tolist())
        dict_frequency_b = dict(frequency_b)

        gene_degree = {}
        for i in gene_set:
            degree_a = 0 if dict_frequency_a.get(i) is None else dict_frequency_a.get(i)
            degree_b = 0 if dict_frequency_b.get(i) is None else dict_frequency_b.get(i)
            gene_degree[i] = degree_a + degree_b
        df_gene_degree_all = pd.DataFrame(gene_degree.items(), columns=['gene','degree'])
        df_gene_degree_all.sort_values(by=['degree'], inplace=True, ascending=False)
        df_gene_degree_all = df_gene_degree_all.set_index('gene')
        
        # add degree score
        df_disease_phewas_or_gwas = pd.merge(df_disease_phewas_or_gwas, df_gene_degree_all, left_on='gene_name', right_on = "gene", how='left')

        # generate PPI network
        for i in df_protein_interaction.index:
            g_all.add_node(df_protein_interaction["preferredName_A"][i], color = "rgb(203,213,232)")
            g_all.add_node(df_protein_interaction["preferredName_B"][i], color = "rgb(203,213,232)")
            g_all.add_edge(df_protein_interaction["preferredName_A"][i],df_protein_interaction["preferredName_B"][i], weight= df_protein_interaction["score"][i])
    except Exception:
        pass

    try:
        st.header("Disease: " + "*" + select_disease + "*" + ", P-value <= " + "*" + str(round(select_p, 4)) + "*")

        with st.container():
            st.subheader("Overall protein-protein interaction network")
            col1, col2, col3= st.columns([5, 1.5, 3.5])
            with col1:
                pv_static(g_all)
            with col2:
                st.write("Protein network degree")
                st.dataframe(df_gene_degree_all, height = 520)
            with col3:
                st.write("GO enrichment analysis")
                st.dataframe(go_enrichment(input_gene = df_disease_phewas_or_gwas['gene_name'].tolist()), height = 520)

        with st.container():
            st.subheader("Risk allele protein-protein interaction network (odds ratio > 1)")
            col1, col2, col3= st.columns([5, 1.5, 3.5])
            with col1:
                pv_static(g_des)
            with col2:
                st.write("Protein network degree")
                st.dataframe(df_gene_degree_des, height = 520)
            with col3:
                st.write("GO enrichment analysis")
                st.dataframe(go_enrichment(input_gene = df_disease_phewas_or_gwas_des['gene_name'].tolist()), height = 520)

        with st.container():
            st.subheader("Protective allele protein-protein interaction network (odds ratio < 1)")
            col1, col2, col3= st.columns([5, 1.5, 3.5])
            with col1:
                pv_static(g_pro)
            with col2:
                st.write("Protein network degree")
                st.dataframe(df_gene_degree_pro, height = 520)
            with col3:
                st.write("GO enrichment analysis")
                st.dataframe(go_enrichment(input_gene = df_disease_phewas_or_gwas_pro['gene_name'].tolist()), height = 520)
    except:
        st.subheader("No data found. Please try selecting another disease!")

elif select == "Disease Target Prioritization":
    st.markdown("This dashboard shows the overall score for target prioritisation of genes with associations with your disease of interest. Features that contribute to the overall score include detection of association and association odds-ratio from a variety of studies (**OpenTargets** and the **PheWAS catalog**), network degree in protein-protein interaction networks (**STRING**), and level of druggability (**Pharos**). ")

    try:
        # extract diseases
        disease = ["--"]
        for i in list(set(df_selected["gwas-associations"].tolist())):
            disease.extend(i.split(", "))
        disease.extend(df_selected["phewas phenotype"].tolist())
        disease = sorted(list(set(disease)))

        # sidebar -- disease select box
        select_disease = st.sidebar.selectbox('Disease', disease, key='11')
        select_disease = select_disease.split(" (")[0]

        # sidebar --  p-value slider
        select_p = st.sidebar.text_input(label = "P-value", help = "Defaults to p = 0.05. Accepts scientific notation, e.g., 5E-4, 3e-9", value = "0.05")
        try:
            if (float(select_p) <= 1) & (float(select_p) > 0):
                select_p = float(select_p)
            else:
                select_p = 0.05
        except:
            select_p = 0.05
        df_selected = df_selected[df_selected["p-value"] <= select_p]
        
        # opentargets gene - disease association score
        try:
            df_gene_score = opentargets_gene_score(disease_name= select_disease)
            #st.write(df_gene_score)

            # add the association score to phewas dataset
            df_selected = pd.merge(df_selected, df_gene_score, left_on='gene_name', right_on = "gene_symbol", how='left')
            df_selected['opentargets_associations'] = df_selected['opentargets_associations'].fillna(0)
        except Exception:
            df_selected['opentargets_associations'] = 0
        #st.write(df_selected['opentargets_associations'])

        # subset the data frame for the disease appearing in either PheWAS and GWAS
        df_disease_phewas_or_gwas = df_selected[df_selected["phewas phenotype"].str.contains(select_disease) | df_selected["gwas-associations"].str.contains(select_disease)]
        df_disease_phewas_or_gwas.sort_values(by=['odds-ratio'], inplace=True, ascending=False)
        df_disease_phewas_or_gwas = df_disease_phewas_or_gwas[["gene_name", "snp", "phewas phenotype", "odds-ratio", "p-value", "gwas-associations", "opentargets_associations"]].reset_index().drop("index", axis= 1)

        # add phewas and gwas interaction indicator score look for the genes in the intersection of phewas and gwas association
        df_disease_phewas_or_gwas = df_disease_phewas_or_gwas.assign(indicator_Phe_GWAS = np.where(df_disease_phewas_or_gwas["phewas phenotype"].str.contains(select_disease) & df_disease_phewas_or_gwas["gwas-associations"].str.contains(select_disease), 1, 0))

        # look for evidence of druggability
        # druggable evidence
        df_disease_phewas_or_gwas = pd.merge(df_disease_phewas_or_gwas, df_druggable, left_on='gene_name', right_on = "sym", how='left')
        
        # subset the data by odds ratio
        df_disease_phewas_or_gwas = df_disease_phewas_or_gwas.reset_index().drop("index", axis= 1)
        df_disease_phewas_or_gwas_des = df_disease_phewas_or_gwas[df_disease_phewas_or_gwas["odds-ratio"] >= 1]
        df_disease_phewas_or_gwas_pro = df_disease_phewas_or_gwas[df_disease_phewas_or_gwas["odds-ratio"] < 1]

        # add druggable score
        # risk gene
        df_tdl_des = df_disease_phewas_or_gwas_des[["gene_name", "tdl"]]
        df_tdl_des =df_tdl_des.dropna()
        df_tdl_des = df_tdl_des.drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
        df_tdl_des = df_tdl_des[df_tdl_des["tdl"] != 'Tdark']
        df_tdl_des_degree = df_tdl_des.groupby(["gene_name"]).size().reset_index(name='druggability_score')
        df_disease_phewas_or_gwas_des = pd.merge(df_disease_phewas_or_gwas_des, df_tdl_des_degree, left_on='gene_name', right_on = "gene_name", how='left')
        df_disease_phewas_or_gwas_des["druggability_score"] = df_disease_phewas_or_gwas_des["druggability_score"].fillna(0)
        
        # protective gene
        df_tdl_pro = df_disease_phewas_or_gwas_pro[["gene_name", "tdl"]]
        df_tdl_pro =df_tdl_pro.dropna()
        df_tdl_pro = df_tdl_pro.drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
        df_tdl_pro = df_tdl_pro[df_tdl_pro["tdl"] != 'Tdark']
        df_tdl_pro_degree = df_tdl_pro.groupby(["gene_name"]).size().reset_index(name='druggability_score')
        df_disease_phewas_or_gwas_pro = pd.merge(df_disease_phewas_or_gwas_pro, df_tdl_pro_degree, left_on='gene_name', right_on = "gene_name", how='left')
        df_disease_phewas_or_gwas_pro["druggability_score"] = df_disease_phewas_or_gwas_pro["druggability_score"].fillna(0)

        # overall gene
        df_tdl_all = df_disease_phewas_or_gwas[["gene_name", "tdl"]]
        df_tdl_all =df_tdl_all.dropna()
        df_tdl_all = df_tdl_all.drop_duplicates(keep = 'first').reset_index().drop('index', axis=1)
        df_tdl_all = df_tdl_all[df_tdl_all["tdl"] != 'Tdark']
        df_tdl_all_degree = df_tdl_all.groupby(["gene_name"]).size().reset_index(name='druggability_score')
        df_disease_phewas_or_gwas = pd.merge(df_disease_phewas_or_gwas, df_tdl_all_degree, left_on='gene_name', right_on = "gene_name", how='left')
        df_disease_phewas_or_gwas["druggability_score"] = df_disease_phewas_or_gwas["druggability_score"].fillna(0)

        # protein-protein network (des genes)
        gene_set = df_disease_phewas_or_gwas_des["gene_name"].unique().tolist()
        g_des=net.Network()
        try:
            df_protein_interaction = proteins_interaction(input_protein= gene_set)
            
            # compute the degree of nodes
            frequency_a = collections.Counter(df_protein_interaction["preferredName_A"].tolist())
            dict_frequency_a = dict(frequency_a)
            frequency_b = collections.Counter(df_protein_interaction["preferredName_B"].tolist())
            dict_frequency_b = dict(frequency_b)

            gene_degree = {}
            for i in gene_set:
                degree_a = 0 if dict_frequency_a.get(i) is None else dict_frequency_a.get(i)
                degree_b = 0 if dict_frequency_b.get(i) is None else dict_frequency_b.get(i)
                gene_degree[i] = degree_a + degree_b
            df_gene_degree = pd.DataFrame(gene_degree.items(), columns=['gene', 'networkDegree_score'])
            
            # add degree score
            df_disease_phewas_or_gwas_des = pd.merge(df_disease_phewas_or_gwas_des, df_gene_degree, left_on='gene_name', right_on = "gene", how='left')

            # generate PPI network
            for i in df_protein_interaction.index:
                g_des.add_node(df_protein_interaction["preferredName_A"][i], color = "#FF7F7F")
                g_des.add_node(df_protein_interaction["preferredName_B"][i], color = "#FF7F7F")
                g_des.add_edge(df_protein_interaction["preferredName_A"][i],df_protein_interaction["preferredName_B"][i], weight= df_protein_interaction["score"][i])
        except Exception:
            df_disease_phewas_or_gwas_des['networkDegree_score'] = 0
            pass

        # protein-protein network (pro genes)
        gene_set = df_disease_phewas_or_gwas_pro["gene_name"].unique().tolist()
        g_pro=net.Network()
        try:
            df_protein_interaction = proteins_interaction(input_protein= gene_set)
            
            # compute the degree of nodes
            frequency_a = collections.Counter(df_protein_interaction["preferredName_A"].tolist())
            dict_frequency_a = dict(frequency_a)
            frequency_b = collections.Counter(df_protein_interaction["preferredName_B"].tolist())
            dict_frequency_b = dict(frequency_b)

            gene_degree = {}
            for i in gene_set:
                degree_a = 0 if dict_frequency_a.get(i) is None else dict_frequency_a.get(i)
                degree_b = 0 if dict_frequency_b.get(i) is None else dict_frequency_b.get(i)
                gene_degree[i] = degree_a + degree_b
            df_gene_degree = pd.DataFrame(gene_degree.items(), columns=['gene', 'networkDegree_score'])
            
            # add degree score
            df_disease_phewas_or_gwas_pro = pd.merge(df_disease_phewas_or_gwas_pro, df_gene_degree, left_on='gene_name', right_on = "gene", how='left')

            # generate PPI network
            for i in df_protein_interaction.index:
                g_pro.add_node(df_protein_interaction["preferredName_A"][i], color = "#45b6fe")
                g_pro.add_node(df_protein_interaction["preferredName_B"][i], color = "#45b6fe")
                g_pro.add_edge(df_protein_interaction["preferredName_A"][i],df_protein_interaction["preferredName_B"][i], weight= df_protein_interaction["score"][i])
        except Exception:
            df_disease_phewas_or_gwas_pro['networkDegree_score'] = 0
            pass

        # protein-protein network (overall genes)
        gene_set = df_disease_phewas_or_gwas["gene_name"].unique().tolist()
        g_all=net.Network()
        try:
            df_protein_interaction = proteins_interaction(input_protein= gene_set)
            
            # compute the degree of nodes
            frequency_a = collections.Counter(df_protein_interaction["preferredName_A"].tolist())
            dict_frequency_a = dict(frequency_a)
            frequency_b = collections.Counter(df_protein_interaction["preferredName_B"].tolist())
            dict_frequency_b = dict(frequency_b)

            gene_degree = {}
            for i in gene_set:
                degree_a = 0 if dict_frequency_a.get(i) is None else dict_frequency_a.get(i)
                degree_b = 0 if dict_frequency_b.get(i) is None else dict_frequency_b.get(i)
                gene_degree[i] = degree_a + degree_b
            df_gene_degree = pd.DataFrame(gene_degree.items(), columns=['gene', 'networkDegree_score'])
            
            # add degree score
            df_disease_phewas_or_gwas = pd.merge(df_disease_phewas_or_gwas, df_gene_degree, left_on='gene_name', right_on = "gene", how='left')

            # generate PPI network
            for i in df_protein_interaction.index:
                g_all.add_node(df_protein_interaction["preferredName_A"][i], color = "#45b6fe")
                g_all.add_node(df_protein_interaction["preferredName_B"][i], color = "#45b6fe")
                g_all.add_edge(df_protein_interaction["preferredName_A"][i],df_protein_interaction["preferredName_B"][i], weight= df_protein_interaction["score"][i])
        except Exception:
            df_disease_phewas_or_gwas['networkDegree_score'] = 0
            pass
        
        # option to remove features from framework
        featuresList = ["odds-ratio", "opentargets_associations", "indicator_Phe_GWAS", "druggability_score", "networkDegree_score"]
        feature_remove = st.sidebar.multiselect("Remove: ", featuresList, key='2')
        if len(feature_remove) > 0:
            featuresList = [ele for ele in featuresList if ele not in feature_remove]
        
        weight = 1/len(featuresList)
        # group the data by gene and normalize features (risk)
        df_disease_phewas_or_gwas_des = df_disease_phewas_or_gwas_des[["gene_name", "odds-ratio", "opentargets_associations", "indicator_Phe_GWAS", "druggability_score", "networkDegree_score"]]
        df_disease_phewas_or_gwas_des = df_disease_phewas_or_gwas_des.groupby(["gene_name"]).mean()
        # normalize by features
        scaler = preprocessing.MinMaxScaler()
        names = df_disease_phewas_or_gwas_des.columns
        ind = df_disease_phewas_or_gwas_des.index
        d_des = scaler.fit_transform(df_disease_phewas_or_gwas_des)
        df_disease_phewas_or_gwas_des_norm = pd.DataFrame(d_des, columns = names, index= ind)
        #df_disease_phewas_or_gwas_des_norm["overall score"] = weight * (df_disease_phewas_or_gwas_des_norm["odds-ratio"] + df_disease_phewas_or_gwas_des_norm["opentargets_associations"] + df_disease_phewas_or_gwas_des_norm["indicator_Phe_GWAS"] + df_disease_phewas_or_gwas_des_norm["druggability_score"] + df_disease_phewas_or_gwas_des_norm["networkDegree_score"])
        featuresList_des = list(featuresList)
        df_disease_phewas_or_gwas_des_norm = df_disease_phewas_or_gwas_des_norm[featuresList_des]
        df_disease_phewas_or_gwas_des_norm["overall score"] = weight * df_disease_phewas_or_gwas_des_norm.sum(axis = 1, skipna = True)
        #df_disease_phewas_or_gwas_des_norm = df_disease_phewas_or_gwas_des_norm[["overall score", "odds-ratio", "opentargets_associations", "indicator_Phe_GWAS", "druggability_score", "networkDegree_score"]]
        featuresList_des.insert(0,"overall score")
        df_disease_phewas_or_gwas_des_norm = df_disease_phewas_or_gwas_des_norm[featuresList_des]
        df_disease_phewas_or_gwas_des_norm.sort_values(by=['overall score'], inplace=True, ascending=False)

        # group the data by gene and normalize features (protective)
        df_disease_phewas_or_gwas_pro["1-odds_ratio"] = 1- df_disease_phewas_or_gwas_pro["odds-ratio"]
        featuresList_pro =['1-odds_ratio' if i =='odds-ratio' else i for i in featuresList]
        df_disease_phewas_or_gwas_pro = df_disease_phewas_or_gwas_pro[["gene_name", "1-odds_ratio", "opentargets_associations", "indicator_Phe_GWAS", "druggability_score", "networkDegree_score"]]
        df_disease_phewas_or_gwas_pro = df_disease_phewas_or_gwas_pro.groupby(["gene_name"]).mean()
        # normalize by features
        scaler = preprocessing.MinMaxScaler()
        names = df_disease_phewas_or_gwas_pro.columns
        ind = df_disease_phewas_or_gwas_pro.index
        d_pro = scaler.fit_transform(df_disease_phewas_or_gwas_pro)
        df_disease_phewas_or_gwas_pro_norm = pd.DataFrame(d_pro, columns = names, index= ind)
        #df_disease_phewas_or_gwas_pro_norm["overall score"] = weight * (df_disease_phewas_or_gwas_pro_norm["1-odds_ratio"] + df_disease_phewas_or_gwas_pro_norm["opentargets_associations"] + df_disease_phewas_or_gwas_pro_norm["indicator_Phe_GWAS"] + df_disease_phewas_or_gwas_pro_norm["druggability_score"] + df_disease_phewas_or_gwas_pro_norm["networkDegree_score"])
        df_disease_phewas_or_gwas_pro_norm = df_disease_phewas_or_gwas_pro_norm[featuresList_pro]
        df_disease_phewas_or_gwas_pro_norm["overall score"] = weight * df_disease_phewas_or_gwas_pro_norm.sum(axis = 1, skipna = True)
        #df_disease_phewas_or_gwas_pro_norm = df_disease_phewas_or_gwas_pro_norm[["overall score", "1-odds_ratio", "opentargets_associations", "indicator_Phe_GWAS", "druggability_score", "networkDegree_score"]]
        featuresList_pro.insert(0,"overall score")
        df_disease_phewas_or_gwas_pro_norm = df_disease_phewas_or_gwas_pro_norm[featuresList_pro]
        df_disease_phewas_or_gwas_pro_norm.sort_values(by=['overall score'], inplace=True, ascending=False)
        df_disease_phewas_or_gwas_pro_norm = df_disease_phewas_or_gwas_pro_norm.rename(columns = {"overall score": "StarGazer score"})

        
        # group the data by gene and normalize features (overall)
        df_disease_phewas_or_gwas = df_disease_phewas_or_gwas[["gene_name", "odds-ratio", "opentargets_associations", "indicator_Phe_GWAS", "druggability_score", "networkDegree_score"]]
        df_disease_phewas_or_gwas = df_disease_phewas_or_gwas.groupby(["gene_name"]).mean()
        # normalize by features
        scaler = preprocessing.MinMaxScaler()
        names = df_disease_phewas_or_gwas.columns
        ind = df_disease_phewas_or_gwas.index
        d_des = scaler.fit_transform(df_disease_phewas_or_gwas)
        df_disease_phewas_or_gwas_norm = pd.DataFrame(d_des, columns = names, index= ind)
        #df_disease_phewas_or_gwas_norm["overall score"] = weight * (df_disease_phewas_or_gwas_norm["odds-ratio"] + df_disease_phewas_or_gwas_norm["opentargets_associations"] + df_disease_phewas_or_gwas_norm["indicator_Phe_GWAS"] + df_disease_phewas_or_gwas_norm["druggability_score"] + df_disease_phewas_or_gwas_norm["networkDegree_score"])
        featuresList_all = list(featuresList)
        df_disease_phewas_or_gwas_norm = df_disease_phewas_or_gwas_norm[featuresList_all]
        df_disease_phewas_or_gwas_norm["overall score"] = weight * df_disease_phewas_or_gwas_norm.sum(axis = 1, skipna = True)
        #df_disease_phewas_or_gwas_norm = df_disease_phewas_or_gwas_norm[["overall score", "odds-ratio", "opentargets_associations", "indicator_Phe_GWAS", "druggability_score", "networkDegree_score"]]
        featuresList_all.insert(0,"overall score")
        df_disease_phewas_or_gwas_norm = df_disease_phewas_or_gwas_norm[featuresList_all]
        df_disease_phewas_or_gwas_norm.sort_values(by=['overall score'], inplace=True, ascending=False)

        df_disease_phewas_or_gwas_norm = df_disease_phewas_or_gwas_norm.rename(columns = {"overall score": "StarGazer score"})
        st.header("Disease: " + "*" + select_disease + "*" + ", P-value <= " + "*" + str(round(select_p, 4)) + "*")

        # target prioritization data frame
        st.subheader("Overall target prioritization: " + "*" + str(len(df_disease_phewas_or_gwas_norm)) + "* genes")
        with st.container():
            st.markdown(get_table_download_link(df_disease_phewas_or_gwas_norm.reset_index()), unsafe_allow_html=True)
            st.dataframe(df_disease_phewas_or_gwas_norm, width = 1100)

        st.subheader("Risk allele target prioritization: " + "*" + str(len(df_disease_phewas_or_gwas_des_norm)) + "* genes")
        with st.container():
            st.markdown(get_table_download_link(df_disease_phewas_or_gwas_des_norm.reset_index()), unsafe_allow_html=True)
            st.dataframe(df_disease_phewas_or_gwas_des_norm, width = 1100)

        st.subheader("Protective allele target prioritization: " + "*" + str(len(df_disease_phewas_or_gwas_pro_norm)) + "* genes")
        with st.container():
            st.markdown(get_table_download_link(df_disease_phewas_or_gwas_pro_norm.reset_index()), unsafe_allow_html=True)
            st.dataframe(df_disease_phewas_or_gwas_pro_norm, width = 1100)

    except:
        st.subheader("No data found. Please try selecting another disease!")

    cache = ""
    placeholder = st.sidebar.empty()
    input = placeholder.text_input("Search NCBI for your gene of interest")
    if input != cache:
        gene_of_interest_URL = "https://www.ncbi.nlm.nih.gov/gene/?term=" + input + "+AND+human[orgn]"
        webbrowser.open(gene_of_interest_URL)
        cache = input
        

else:
    st.markdown("StarGazer is a multi-omics pipeline which integrates several datasets to provide insights into therapeutic target prioritisation. We have integrated data from [OpenTargets](https://www.opentargets.org/) and the [PheWAS catalog](https://phewascatalog.org/phewas) (gene variant risk associations with phenotypic variants), [Pharos](https://pharos.nih.gov/) (druggability of gene target), and [STRING](https://string-db.org/) (protein-protein interaction).")

    # extract GWAS diseases
    disease = []
    for i in list(set(df_selected["gwas-associations"].tolist())):
        disease.extend(i.split(", "))
    disease = sorted(list(set(disease)))
    num_gwas = len(disease)

    # extract phewas diseases
    disease = []
    disease.extend(df_selected["phewas phenotype"].tolist())
    disease = sorted(list(set(disease)))
    num_phewas = len(disease)

    disease = []
    for i in list(set(df_selected["gwas-associations"].tolist())):
        disease.extend(i.split(", "))
    disease.extend(df_selected["phewas phenotype"].tolist())
    disease = sorted(list(set(disease)))
    num_phewas_gwas = len(disease)

    st.markdown("#### StarGazer provides various functionalities - users can search by:")
    st.markdown("- Gene")
    st.markdown("- Variant")
    st.markdown("- PheWAS")
    st.markdown("- GWAS")
    st.markdown("- GWAS_PheWAS Union")
    st.markdown("- GWAS_PheWAS Intersection")
    st.markdown("- Protein-protein Interaction")
    st.markdown("- Disease Target Prioritization")
    st.markdown("#### Please start your navigation from the sidebar!")
    st.markdown("Users should note that drop-down menus also function as search boxes.")

