# -------------------------------------------------------
# Cleaning data extracted from epigraph db
# -------------------------------------------------------

# This script aims to format the data extracted from epigraphdb
# so that it matches the GO and KEGG term data 

pkgs <- c("tidyverse", "FDb.InfiniumMethylation.hg19",
		  "IlluminaHumanMethylation450kanno.ilmn12.hg19")

lapply(pkgs, require, character.only = T)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

devtools::load_all("~/repos/usefunc/")

source("R/mapping_functions.R")
source("R/simulation_functions.R")

# read in kegg terms for reference
kegg_terms <- read_tsv("data/derived/kegg_terms.txt")

# read in gene data
all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max=1e6)


# read in epigraph db data
reactome_dat <- read_tsv("data/derived/epigraphdb_proteins_and_pathways.tsv")
epidb_proteins <- read_tsv("data/derived/epigraphdb_uniprot_ids.tsv")
epidb_ppi <- read_tsv("data/derived/epigraphdb_proteins-protein_interactions.tsv")

# -------------------------------------------------------
# Clean the data
# -------------------------------------------------------

sum(epidb_proteins$ensembl_gene_id %in% all_genes$ensembl_gene_id)

react_pathway_names <- reactome_dat %>%
	pull(pathway_reactome_id) %>%
	gsub("\\[", "", .) %>%
	gsub("\\]", "", .) %>%
	gsub("\\'", "", .) %>%
	strsplit(", ") %>%
	unlist() %>%
	unique

x=react_pathway_names[1]
reactome_terms <- map_dfr(react_pathway_names, function(x) {
	proteins_in_pathway <- reactome_dat %>%
		dplyr::filter(grepl(x, pathway_reactome_id)) %>%
		pull(uniprot_id) %>%
		unique
	out <- epidb_proteins %>%
		dplyr::filter(uniprot_id %in% proteins_in_pathway) %>%
		mutate(pathway_id = x) %>%
		dplyr::select(pathway_id, ensembl_gene_id)
	return(out)
})

ppi_terms <- epidb_ppi %>%
	dplyr::filter(!is.na(assoc_protein)) %>%
	left_join(epidb_proteins, by = c("assoc_protein" = "uniprot_id")) %>%
	dplyr::select(pathway_id = protein, ensembl_gene_id) %>%
	distinct()

# -------------------------------------------------------
# Write out the data
# -------------------------------------------------------

write.table(reactome_terms, "data/derived/reactome_terms.txt", 
			col.names = T, row.names = F, quote = F, sep = "\t")

write.table(ppi_terms, "data/derived/ppi_terms.txt", 
			col.names = T, row.names = F, quote = F, sep = "\t")