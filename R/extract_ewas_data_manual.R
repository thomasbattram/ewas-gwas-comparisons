# ---------------------------------------------------------
# Script to sort the ewas data that are external to ewas catalog
# ---------------------------------------------------------

## This script links to extract_ewas_data.R and is used to manually
## extract any data that may not be in The EWAS Catalog and can be 
## included in the study. 

## pkgs
library(tidyverse) # tidy code and data

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

devtools::load_all("~/repos/usefunc/")

## read in ewas catalog data for reference
ec_studies <- read_tsv("data/ewas_catalog/studies.txt")
ec_results <- read_tsv("data/ewas_catalog/results.txt", n_max=10)

## read in the external ewas data
glu_ewas <- read_tsv("data/glucose_ewas/FG.model1.txt") # model 1 is without adjustment for BMI
ins_ewas <- read_tsv("data/glucose_ewas/FI.model1.txt") # model 1 is without adjustment for BMI

# ---------------------------------------------------------
# Sort functions for glucose + insulin ewas
# ---------------------------------------------------------

make_studies_template <- function(ec_studies)
{
	### function to make studies template for external data
	studies_template <- matrix(ncol = length(colnames(ec_studies)), nrow = 1) %>%
		as_tibble()
	colnames(studies_template) <- colnames(ec_studies)
	return(studies_template)
}

make_results_template <- function(ec_results, rows)
{
	### function to make studies template for external data
	results_template <- matrix(ncol = length(colnames(ec_results)), nrow = rows) %>%
		as_tibble()
	colnames(results_template) <- colnames(ec_results)
	return(results_template)
}


sort_glucose_studies <- function(ewas_dat, ec_studies, trait)
{
	### function to sort studies for glucose + insulin ewas
	out <- make_studies_template(ec_studies) %>%
		mutate(N = max(ewas_dat$N), 
			   Author = "Liu J", 
			   Trait = trait, 
			   Exposure = trait, 
			   Outcome = "DNA methylation", 
			   PMID = 31197173, 
			   StudyID = paste0("external_ewas_", trait))
}

# var_of_interest <- c("cpg", "beta", "se", "pval", "N")

sort_glucose_results <- function(ewas_dat, ec_results, trait)
{
	### function to sort results for glucose + insulin ewas
	results_template <- make_results_template(ec_results, rows=nrow(ewas_dat))
	out <- make_results_template(ec_results, rows=nrow(ewas_dat)) %>%
		mutate(CpG = ewas_dat$cpg, 
			   Beta = ewas_dat$beta, 
			   SE = ewas_dat$se, 
			   P = ewas_dat$pval, 
			   StudyID = paste0("external_ewas_", trait))
}

# ---------------------------------------------------------
# Sorting glucose + insulin EWAS
# ---------------------------------------------------------

gi_ewas <- list(glucose = glu_ewas, insulin = ins_ewas)
lapply(seq_along(gi_ewas), function(x) {
	trait <- names(gi_ewas)[x]
	ewas_dat <- gi_ewas[[x]]
	studies <- sort_glucose_studies(ewas_dat, ec_studies, trait)
	results <- sort_glucose_results(ewas_dat, ec_results, trait)
	out <- list(studies = studies, results = results)
	out_nam <- paste0("data/glucose_ewas/", trait, "_ewas_data.RData")
	save(out, file = out_nam)
})
