# ---------------------------------------------------------
# Script to sort the ewas data that are external to ewas catalog
# ---------------------------------------------------------

## This script links to extract_ewas_data.R and is used to manually
## extract any data that may not be in The EWAS Catalog and can be 
## included in the study. 

## pkgs
library(tidyverse) # tidy code and data
library(usefunc) # own package of useful functions

## read in ewas catalog data for reference
ec_studies <- read_tsv("data/ewas_catalog/studies.txt", guess_max = 1e6)
ec_results <- read_tsv("data/ewas_catalog/results.txt", n_max=10)

## read in the external ewas data
get_ewas_from_web <- function(out_path, web_path)
{
	if (file.exists(out_path)) {
		out <- read_csv(out_path)
	} else {
		out_dir <- dirname(out_path)
		if (!file.exists(out_dir)) {
			message("Creating directory:")
			message(getwd(), "/", out_dir)
			system(paste0("mkdir -p ", out_dir))
		}
		d_status <- download.file(url = web_path, 
								  destfile = out_path,
								  method = "wget")
		if (d_status != 0) stop("Download did not complete. Please retry or check filepaths.")
		out <- read_csv(out_path)
	}
	return(out)
}

egfr_ewas <- get_ewas_from_web(out_path = "data/ckd_ewas/eGFR.csv.zip", 
						   	   web_path = "https://ckdgen.imbi.uni-freiburg.de/files/Schlosser2021/eGFR.csv.zip")
urate_ewas <- get_ewas_from_web(out_path = "data/ckd_ewas/urate.csv.zip", 
								web_path = "https://ckdgen.imbi.uni-freiburg.de/files/Tin2021/urate.csv.zip")
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
	make_studies_template(ec_studies) %>%
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
	return(out)
}

sort_ckd_studies <- function(ewas_dat, ec_studies, trait)
{
	### function to sort studies for egfr and urate EWAS
	make_studies_template(ec_studies) %>%
		mutate(N = max(ewas_dat$n), 
			   Author = ifelse(trait == "egfr", "Schlosser P", "Tin A"), 
			   Trait = trait, 
			   Exposure = trait, 
			   Outcome = "DNA methylation", 
			   PMID = ifelse(trait == "egfr", 34887417, 34887389), 
			   StudyID = paste0("external_ewas_", trait))
}

sort_ckd_results <- function(ewas_dat, ec_results, trait)
{
	### function to sort results for egfr + urate ewas
	results_template <- make_results_template(ec_results, rows=nrow(ewas_dat))
	out <- make_results_template(ec_results, rows=nrow(ewas_dat)) %>%
		mutate(CpG = ewas_dat$probeID, 
			   Beta = ewas_dat[["Effect estimate"]], 
			   SE = ewas_dat[["Standard error"]], 
			   P = ewas_dat[["P-value"]], 
			   StudyID = paste0("external_ewas_", trait))
	return(out)
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
rm(gi_ewas)

ckd_ewas <- list(egfr = egfr_ewas, urate = urate_ewas)
lapply(seq_along(ckd_ewas), function(x) {
	trait <- names(ckd_ewas)[x]
	ewas_dat <- ckd_ewas[[x]]
	studies <- sort_ckd_studies(ewas_dat, ec_studies, trait)
	results <- sort_ckd_results(ewas_dat, ec_results, trait)
	out <- list(studies = studies, results = results)
	out_nam <- paste0("data/ckd_ewas/", trait, "_ewas_data.RData")
	save(out, file = out_nam)
})

