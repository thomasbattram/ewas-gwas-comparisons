# -------------------------------------------------------
# Examining enriched pathways
# -------------------------------------------------------

## goal of script:
## to assess which pathways GWAS and EWAS genes may be enriched in
## and compare the actual pathway descriptions 

# 1. get empirical genes 
# 2. perform enrichment analyses with these genes
# 3. compare pathways enriched for amongst corresponding GWAS and EWAS

# pkgs
library(tidyverse) # tidy code and data 
library(usefunc) # own package of useful functions

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

source("R/mapping_functions.R")
source("R/simulation_functions.R")

# -------------------------------------------------------
# Read in necessary data
# -------------------------------------------------------

pathway_dat_nam <- "data/derived/go_terms.txt"
pathway_dat <- read_tsv(pathway_dat_nam)

traits <- readLines("data/traits.txt")
cpg_genes <- read_tsv("data/derived/epigenetic/cpgs_mapped_to_genes.txt")
snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_genes.txt")
all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)  %>%
    dplyr::filter(!duplicated(ensembl_gene_id))

path <- "data/derived"
db_gene_file <- "go_genes.rdata"
full_path <- file.path(path, db_gene_file)
if (file.exists(full_path)) {
    load(full_path)
    assign("go_genes", pathway_genes)
} else {
    pathway_df <- get("go_terms")
    pathway_genes <- extract_genes_in_pathways(pathway_df)
    save(pathway_genes, file = full_path)
    assign("go_genes", pathway_genes)
}

# -------------------------------------------------------
# Perform enrichment for each trait
# -------------------------------------------------------

extract_ewas_and_gwas_genes <- function(trait, cpg_genes, snp_genes)
{
	### extracting "significant" genes from ewas and gwas
	
	# read in data
	dat <- read_derived_data(trait)
	# extract ewas genes
	sig_cpgs <- dat$meth[dat$meth$P < 1e-7, "CpG", drop = T]
	ewasg <- unique(cpg_genes[cpg_genes$name %in% sig_cpgs, "ensembl_gene_id", drop = T])
	# extract gwas genes
	sig_snps <- dat$gen[dat$gen$p < 5e-8, "SNP", drop = T]
	gwasg <- unique(snp_genes[snp_genes$name %in% sig_snps, "ensembl_gene_id", drop = T])
	return(list(ewasg = ewasg, gwasg = gwasg))
}

get_pathway_size <- function(pathways, pathway_dat)
{
	pathway_dat %>%
		dplyr::filter(pathway_id %in% pathways) %>%
		group_by(pathway_id) %>%
		dplyr::summarise(genes_in_pathway = n())
}

get_enriched_paths <- function(enrichment_res, pathway_dat)
{
	enrich_out <- enrichment_res %>%
		mutate(adj_enrich_p = p.adjust(enrich_p, method = "BH")) %>%
		# dplyr::filter(adj_enrich_p < 0.1) %>%
		arrange(enrich_p) %>%
		.[1:100,] %>%
		left_join(pathway_dat) %>% 
		dplyr::select(-ensembl_gene_id) %>%
		distinct()
	pathsizes <- get_pathway_size(enrich_out$pathway_id, pathway_dat)
	enrich_out <- enrich_out %>%
		left_join(pathsizes)
	return(enrich_out)
}

## check pathway sizes
all_pathsizes <- get_pathway_size(unique(pathway_dat$pathway_id), pathway_dat)
summary(all_pathsizes)

## plot to see if there is an obvs cutoff! 
# remove super large pathways first
all_pathsizes <- all_pathsizes %>%
	dplyr::filter(genes_in_pathway < 5000)

pathsize_hist <- ggplot(all_pathsizes, aes(x = genes_in_pathway)) +
	geom_histogram(fill="blue", colour="black", binwidth=100) + 
	theme_bw()
ggsave("results/plots/go_term_sizes.pdf", plot=pathsize_hist)

# looks like less than 100 is a decent start
# also should remove genesets that are just 1 gene...
sum(all_pathsizes$genes_in_pathway < 5)

pathway_out <- lapply(seq_along(traits), function(x) {
	tra <- traits[x]
	sig_g <- extract_ewas_and_gwas_genes(tra, cpg_genes, snp_genes)
	# ewas enrichment
	ewas_enrich <- perform_enrichment(identified_genes = sig_g$ewasg, 
                                  database = "go", 
                                  background_genes = all_genes$ensembl_gene_id)
	ewas_paths <- get_enriched_paths(ewas_enrich, pathway_dat)
	# gwas enrichment
	gwas_enrich <- perform_enrichment(identified_genes = sig_g$gwasg, 
                                  database = "go", 
                                  background_genes = all_genes$ensembl_gene_id)

	gwas_paths <- get_enriched_paths(gwas_enrich, pathway_dat)
	# overlapping paths
	overlapping_path_n <- tibble(n_ewas_paths = nrow(ewas_paths), 
								 n_gwas_paths = nrow(gwas_paths), 
								 n_overlapping_paths = sum(ewas_paths$pathway_id %in% gwas_paths$pathway_id))
	return(list(gwas = gwas_paths, ewas = ewas_paths, summ = overlapping_path_n))
})
names(pathway_out) <- traits

save(pathway_out, file="results/pathway_description_overlap.RData")

# -------------------------------------------------------
# Quick look at enrichment results
# -------------------------------------------------------

pathway_out <- new_load("results/pathway_description_overlap.RData")

pathway_out %>%
	map("summ")

pathway_out %>%
	map("gwas") %>%
	map(function(x) {dplyr::filter(x, adj_enrich_p < 0.1, genes_in_pathway < 100)}) 

pathway_out %>%
	map("gwas") %>%
	map(function(x) {dplyr::filter(x, adj_enrich_p < 0.1, genes_in_pathway < 100)}) %>%
	bind_rows(.id = "trait") %>%
	as.data.frame

pathway_out %>%
	map("ewas") %>%
	map(function(x) {dplyr::filter(x, adj_enrich_p < 0.1, genes_in_pathway < 100)})

pathway_out %>%
	map("ewas") %>%
	map(function(x) {dplyr::filter(x, adj_enrich_p < 0.1, genes_in_pathway < 100)}) %>%
	bind_rows(.id = "trait") %>%
	as.data.frame

# -------------------------------------------------------
# Extract enrichments common (and uncommon) to corresponding traits 
# -------------------------------------------------------

### want to have something in the paper like this:
## There were XXX GO term genesets that were commonly enriched (FDR < 0.1) for both the 
## EWAS and GWAS traits. Of these XXX were non-specific and contained over 100 genes each. There
## were XXX specific genesets (geneset size < 100 genes) that did not overlap between studies of corresponding traits, 
## for example, the genes identified by the GWAS of alcohol consumption were enriched for the 
## “ethanol catabolism” pathway (geneset size = 12 genes), however none of
## the genes identified by the EWAS were present in this pathway. 

## extract pathways common to both traits! 
common_p <- lapply(traits, function(trait) {
	ewas_df <- map(pathway_out, "ewas")[[trait]]
	gwas_df <- map(pathway_out, "gwas")[[trait]]
	all_df <- bind_rows(list(ewas = ewas_df, gwas = gwas_df), .id="study")
	all_df <- dplyr::filter(all_df, adj_enrich_p < 0.1)
	dup_terms <- all_df$pathway_id[duplicated(all_df$pathway_id)]
	if (length(dup_terms) == 0) return(NULL)
	all_df %>%
		dplyr::filter(pathway_id %in% dup_terms)
})
names(common_p) <- traits

all_common <- sum(map_dbl(common_p, function(x) {ifelse(!is.null(x), nrow(x), 0)}))

## now just specific genesets
common_p_specifc <- lapply(common_p, function(df) {
	if (is.null(df)) return(NULL)
	df %>%
		dplyr::filter(genes_in_pathway < 100)
})
specific_common <- sum(map_dbl(common_p_specifc, function(x) {ifelse(!is.null(x), nrow(x), 0)}))
# 0

## extract uncommon pathways 
uncommon_p <- lapply(traits, function(trait) {
	ewas_df <- map(pathway_out, "ewas")[[trait]]
	gwas_df <- map(pathway_out, "gwas")[[trait]]
	all_df <- bind_rows(list(ewas = ewas_df, gwas = gwas_df), .id="study")
	all_df <- dplyr::filter(all_df, adj_enrich_p < 0.1)
	dup_terms <- all_df$pathway_id[duplicated(all_df$pathway_id)]
	all_df %>%
		dplyr::filter(!pathway_id %in% dup_terms)
})
names(uncommon_p) <- traits

all_uncommon <- sum(map_dbl(uncommon_p, function(x) {ifelse(!is.null(x), nrow(x), 0)}))

## now just specific uncommon genesets
uncommon_p_specifc <- lapply(uncommon_p, function(df) {
	if (is.null(df)) return(NULL)
	df %>%
		dplyr::filter(genes_in_pathway < 100)
})
specific_uncommon <- sum(map_dbl(uncommon_p_specifc, function(x) {ifelse(!is.null(x), nrow(x), 0)}))
# 261

## now try and find an intuitive example of differences in "biological"
## understanding gained from different studies

uncommon_p_specifc ## alcohol consumption -- ethanol catabolism

pathway_out$alcohol_consumption_per_day$ewas %>%
	dplyr::filter(adj_enrich_p < 0.1) %>%
	data.frame

pathway_out$alcohol_consumption_per_day$gwas %>%
	dplyr::filter(adj_enrich_p < 0.1) %>%
	data.frame

## quick check for ethanol catabolism in alcohol consumption ewas 

eth_go_id <- "GO:0006068"
	
tra <- "alcohol_consumption_per_day"
sig_g <- extract_ewas_and_gwas_genes(tra, cpg_genes, snp_genes)
# ewas enrichment
ewas_enrich <- perform_enrichment(identified_genes = sig_g$ewasg, 
                            	  database = "go", 
                            	  background_genes = all_genes$ensembl_gene_id)

ewas_enrich %>%
	dplyr::filter(pathway_id == eth_go_id)

# -------------------------------------------------------
# save out things for paper
# -------------------------------------------------------

enriched_pathways_out <- list(n_common = all_common, 
							  n_common_specific = specific_common, 
							  n_uncommon = all_uncommon,
							  n_uncommon_specific = specific_uncommon)

save(enriched_pathways_out, file = "report/report_data/enriched_geneset_examples.RData")
