# ------------------------------------------------------------------------
# Exploring gene overlap
# ------------------------------------------------------------------------

# This script is to assess whether gene overlap tends to be 
# with protein coding genes or not 

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = T)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

devtools::load_all("~/repos/usefunc/")

# source in functions for:
# - reading in the data
# - mapping variant positions to genes
source("R/mapping_functions.R")

#
# test with just BMI data
#

# read in data
cpg_genes <- read_tsv("data/derived/epigenetic/cpgs_mapped_to_genes.txt")
cpg_proteins <- read_tsv("data/derived/epigenetic/cpgs_mapped_to_proteins.txt")
snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_genes.txt")
snp_proteins <- read_tsv("data/derived/genetic/snps_mapped_to_proteins.txt")
all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)
protein_mapping <- read_tsv("data/derived/epigraphdb_uniprot_ids.tsv")

# --------------------------------------------------------
# Extract overlapping genes
# --------------------------------------------------------

extract_sig <- function(derived_data, gen_or_meth, genes_or_proteins)
{
	### Extract sig hits
	stopifnot(gen_or_meth %in% c("gen", "meth"))
	stopifnot(genes_or_proteins %in% c("genes", "proteins"))
	sigp <- ifelse(gen_or_meth =="gen", 5e-8, 1e-7)
	variant_name <- ifelse(gen_or_meth == "gen", "snp", "cpg")
	study_type <- ifelse(gen_or_meth == "gen", "gwas", "ewas")
	df <- derived_data[[gen_or_meth]]
	colnames(df) <- tolower(colnames(df))
	out <- df %>%
		dplyr::filter(p < sigp) %>%
		dplyr::select(name = one_of(variant_name)) %>%
		left_join(get(paste0(variant_name, "_", genes_or_proteins))) %>%
		dplyr::select(ensembl_gene_id, hgnc_symbol) %>%
		mutate(study = study_type)
	return(out)
}

extract_identified_genes <- function(trait)
{
	### For each trait extract the genes that was identified 
	### in the EWAS and GWAS at typical signicance level
	###
	dat <- read_derived_data(trait)
	g_or_p <- c("genes", "proteins")
	ewas_sig <- map_dfr(g_or_p, function(x) {
		out <- extract_sig(dat, gen_or_meth="meth", genes_or_proteins=x)
		out$mapping <- x
		return(out)
	})
	gwas_sig <- map_dfr(g_or_p, function(x) {
		out <- extract_sig(dat, gen_or_meth="gen", genes_or_proteins=x)
		out$mapping <- x
		return(out)
	})	
	out_genes <- rbind(ewas_sig, gwas_sig) %>%
		distinct()
	return(out_genes)
}

traits <- readLines("data/traits.txt")
trait_genes <- lapply(traits, function(tr) {
	message(tr)
	extract_identified_genes(tr)
})
names(trait_genes) <- traits

# --------------------------------------------------------
# Check whether there are many ensembl IDs that don't map to proteins
# --------------------------------------------------------

no_protein_dat <- bind_rows(trait_genes, .id = "trait") %>%
	left_join(protein_mapping) %>%
	dplyr::filter(mapping == "genes") %>%
	group_by(trait, study) %>%
	summarise(prot_na_count = sum(is.na(uniprot_id)), gene_count = n()) %>%
	mutate(proportion_missing = prot_na_count/gene_count)

no_protein_dat
## Fair proportion of genes seem to be mapped to ensembl IDs and NOT protein coding genes

# --------------------------------------------------------
# Check overlap differences between gene and protein mapping
# --------------------------------------------------------

get_overlap <- function(df)
{
	ewas_genes <- df %>%
		dplyr::filter(study == "ewas") %>%
		pull(ensembl_gene_id) %>%
		unique
	gwas_genes <- df %>%
		dplyr::filter(study == "gwas") %>%
		pull(ensembl_gene_id) %>%
		unique
	return(length(intersect(ewas_genes, gwas_genes)))
}

overlap_dat <- map_dfr(traits, function(tr) {
	trait_dat <- trait_genes[[tr]]
	gene_mapping_overlap <- trait_dat %>%
		dplyr::filter(mapping == "genes") %>%
		get_overlap
	protein_mapping_overlap <- trait_dat %>%
		dplyr::filter(mapping == "proteins") %>%
		get_overlap
	out <- tibble(mapping = c("genes", "proteins"), 
				  trait = tr, 
				  overlap = c(gene_mapping_overlap, protein_mapping_overlap))
	return(out)
})

overlap_dat

# --------------------------------------------------------
# Check differences in distance to gene
# --------------------------------------------------------

df_nams <- c("snp_genes", "snp_proteins", 
			 "cpg_genes", "cpg_proteins")

# n within gene
lapply(df_nams, function(x) {
	df <- get(x) 
	out <- data.frame(within_gene_n = sum(df$distance_to_gene == 0),
					  outside_gene_n = sum(df$distance_to_gene > 0))
	out
})

# distance between variant and gene
lapply(df_nams, function(x) {
	df <- get(x) %>%
		dplyr::filter(distance_to_gene > 0)
	summary(df$distance_to_gene)
})
# Median increases by quite a few kb...


