# --------------------------------------------------------
# Extracting epigraph db data using the api
# --------------------------------------------------------

# aim of script: to use epigraphdb api to map genes to proteins 
# and pathways. Should be able to extract  

# What I need at end of it:
# df or lists of genes + linked pathways + protein-protein interactions

from itertools import permutations
from typing import List

# import matplotlib.pyplot as plt
import networkx as nx
import requests, os
import pandas as pd

# setwd
os.getcwd()
os.chdir('main_project/epi_gen_comp')

# connect to api
API_URL = "https://api.epigraphdb.org"
requests.get(f"{API_URL}/ping").json()

# read in genes
all_genes = pd.read_csv("data/derived/ensembl_genes.txt", delimiter="\t")

## map genes to proteins
def get_gene_protein(genelist, ensembl=False):
	""" Map genes to proteins via epigraphdb api
		genes must be hgnc or ensembl gene ids
	""" 
	endpoint = "/mappings/gene-to-protein"
	url = f"{API_URL}{endpoint}"
	if ensembl:
		pl = 'gene_id_list'
	else:
		pl = 'gene_name_list'
	payload = {pl: genelist, 'by_gene_id': ensembl}
	r = requests.post(url, json=payload)
	r.raise_for_status()
	protein_df = pd.json_normalize(r.json()["results"])
	if len(protein_df) > 0:
		res_df = protein_df[["gene.name", "gene.ensembl_id", "protein.uniprot_id"]].rename(
			columns={
				"gene.name": "protein_name",
				"gene.ensembl_id": "ensembl_gene_id",
				"protein.uniprot_id": "uniprot_id",
			}
		)
	else:
		return pd.DataFrame(columns=["protein_name", "uniprot_id"])
	return res_df
 
# GENELIST=list(all_genes['hgnc_symbol'])
GENELIST=list(all_genes['ensembl_gene_id'])
gene_protein_df = get_gene_protein(genelist=GENELIST, ensembl=True)
gene_protein_df

## map proteins to pathways
def get_protein_pathway(gene_protein_df):
	""" Map proteins to reactomedb pathways via epigraphdb api	
	""" 
	endpoint = "/protein/in-pathway"
	url = f"{API_URL}{endpoint}"
	payload = {
		"uniprot_id_list": gene_protein_df["uniprot_id"].to_list(),
	}
	r = requests.post(url, json=payload)
	r.raise_for_status()
	df = pd.json_normalize(r.json()["results"])
	if len(df) > 0:
		res_df = gene_protein_df[["uniprot_id"]].merge(
			df, left_on="uniprot_id", right_on="uniprot_id", how="left"
		)
	else:
		res_df = (
			gene_protein_df[["uniprot_id"]]
			.assign(pathway_count=None)
			.assign(pathway_reactome_id=None)
		)
	res_df = res_df.assign(
		pathway_count=lambda df: df["pathway_count"]
		.apply(lambda x: 0 if pd.isna(x) else x)
		.astype(int)
	).assign(
		pathway_reactome_id=lambda df: df["pathway_reactome_id"].apply(
			lambda x: [] if not isinstance(x, list) else x
		)
	)
	return res_df
 
pathway_df = get_protein_pathway(gene_protein_df=gene_protein_df)
pathway_df

## map protein-protein interactions
def get_ppi(gene_protein_df, n_intermediate_proteins: int = 0):
	""" Map protein-protein interactions via epigraphdb api
		which uses mixture of stringdb and intact
	"""
	endpoint = "/protein/ppi/pairwise"
	url = f"{API_URL}{endpoint}"
	payload = {
		"uniprot_id_list": gene_protein_df["uniprot_id"].to_list(),
		"n_intermediate_proteins": n_intermediate_proteins,
	}
	r = requests.post(url, json=payload)
	r.raise_for_status()
	df = pd.json_normalize(r.json()["results"])
	if len(df) > 0:
		res_df = (
			gene_protein_df[["uniprot_id"]]
			.rename(columns={"uniprot_id": "protein"})
			.merge(df, left_on="protein", right_on="protein", how="left")
		)
	else:
		res_df = (
			gene_protein_df[["uniprot_id"]]
			.rename(columns={"uniprot_id": "protein"})
			.assign(assoc_protein=None, path_size=None)
		)
	return res_df
 
PPI_N_INTERMEDIATE_PROTEINS=0
ppi_df = get_ppi(
	gene_protein_df, n_intermediate_proteins=PPI_N_INTERMEDIATE_PROTEINS
)
ppi_count = ppi_df.groupby(["protein", "assoc_protein"]).count()
ppi_count.rename(columns={"path_size": "count"})
ppi_count
ppi_count.path_size.value_counts()

## output it all! 
out_path = "data/derived/"
# first test pathway df
pathway_df.pathway_count.value_counts()

out_file = "epigraphdb_proteins_and_pathways.tsv"
pathway_df.to_csv(f"{out_path}{out_file}", 
				  index = False, 
				  sep = '\t')

# now output protein-protein interaction data
out_file = "epigraphdb_proteins-protein_interactions.tsv"
ppi_df.to_csv(f"{out_path}{out_file}", 
			  index = False, 
			  sep = '\t')

# finally output proteins linked to uniprot_id
out_file = "epigraphdb_uniprot_ids.tsv"
gene_protein_df.to_csv(f"{out_path}{out_file}", 
			  index = False, 
			  sep = '\t')




