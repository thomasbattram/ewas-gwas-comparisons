# ------------------------------------------------------------
# Mapping CpGs and SNPs to genes and pathways
# ------------------------------------------------------------

## This script aims to take CpGs and SNPs and map them to genes
## based on chromosomal position (nearest gene(s)). Then map
## these genes to genesets

## pkgs
library(tidyverse) # tidy data and code
library(meffil) # 450k probe annotations

args <- commandArgs(trailingOnly = TRUE)
proteins <- as.logical(args[1])

message("Mapping to proteins (T or F): ", proteins)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

devtools::load_all("~/repos/usefunc/")

# source in functions for:
# - reading in the data
# - mapping variant positions to genes
source("R/mapping_functions.R")


# read in gene dat
all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max=1e6)
protein_mapping <- read_tsv("data/derived/epigraphdb_uniprot_ids.tsv")

out_dir <- "data/derived"

# ------------------------------------------------------------
# mapping snp and cpg positions to genes
# ------------------------------------------------------------
if (proteins) {
	gene_mapping_dat <- all_genes %>%
		dplyr::filter(ensembl_gene_id %in% protein_mapping$ensembl_gene_id) %>%
		dplyr::filter(!duplicated(ensembl_gene_id))
	meth_filenam <- file.path(out_dir, "epigenetic", "cpgs_mapped_to_proteins.txt")
	gen_filenam <- file.path(out_dir, "genetic", "snps_mapped_to_proteins.txt")
} else {
	gene_mapping_dat <- all_genes %>%
		dplyr::filter(!duplicated(ensembl_gene_id))
	meth_filenam <- file.path(out_dir, "epigenetic", "cpgs_mapped_to_genes.txt")
	gen_filenam <- file.path(out_dir, "genetic", "snps_mapped_to_genes.txt")
}


## DNA methylation data
mef_probe_info <- meffil.get.features()
head(mef_probe_info)
unique(mef_probe_info$chromosome)
sum(is.na(mef_probe_info$chromosome))
head(mef_probe_info[is.na(mef_probe_info$chromosome), ])

mef_probe_info <- mef_probe_info %>%
	dplyr::filter(grepl("^cg[0-9]*", name)) %>%
	mutate(chromosome = gsub("^chr", "", chromosome)) %>%
	dplyr::select(name, chromosome, position)

# now map the CpG positions to genes

## for tests: mef_probe_info <- mef_probe_info[mef_probe_info$chromosome == "Y", ]

# loop the function over each chromosome
chrom <- unique(mef_probe_info$chromosome)
mapped_dat <- lapply(chrom, function(chr) {
	map_over_chromosomes(chr, variant_dat = mef_probe_info, 
						 gene_dat = gene_mapping_dat)
})
mapped_dat2 <- bind_rows(mapped_dat)

# sanity check -- mapped all of the cpg sites???
length(unique(mapped_dat2$name)) == length(unique(mef_probe_info$name)) # TRUE

write.table(mapped_dat2, file = meth_filenam,
			row.names = F, col.names = T, quote = F, sep = "\t")

# now need to do the same for the genetic data!!! 
# just do it for the divided genome dat data! 
div_dat <- read_tsv("data/derived/genetic/genetic_data_divided_genome.txt")
head(div_dat)
div_dat <- div_dat %>%
	dplyr::rename(position = bp, name = ID) %>%
	dplyr::select(-group, -positions, -group_total) %>%
	dplyr::mutate(CHR = as.character(CHR))

## for tests: div_dat <- div_dat[div_dat$CHR == "22", ]

ori_time <- proc.time()
chrom <- unique(div_dat$CHR)
gen_mapped_dat <- lapply(chrom, function(chr) {
	map_over_chromosomes(chr, variant_dat = div_dat, gene_dat = gene_mapping_dat)
})
gen_mapped_dat <- bind_rows(gen_mapped_dat)
fin_time <- proc.time() - ori_time 
fin_time # 18 mins for 394080 rows
length(unique(gen_mapped_dat$name)) == length(unique(div_dat$name)) # TRUE

write.table(gen_mapped_dat, file = gen_filenam,
			row.names = F, col.names = T, quote = F, sep = "\t")

rm(div_dat)
# FIN



