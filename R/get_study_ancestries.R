# -----------------------------------------------------------------
# Extract ancestries of participants in studies used 
# -----------------------------------------------------------------

## Aim: To extract ancestries from opengwas and ewas catalog 

## pkgs
library(tidyverse) # tidy code and data
library(ieugwasr) # to get gwas data
library(usefunc) # own package of useful functions

## data (assuming in epi_gen_comp directory)
ewas_study_info <- read_tsv("report/report_data/ewas_studies_info.tsv")
gwas_study_info <- read_tsv("report/report_data/gwas_studies_info.tsv")

opengwas_dat <- gwasinfo()
ewascatalog_dat <- read_tsv("http://www.ewascatalog.org/static//docs/ewascatalog-studies.txt.gz", guess_max = 1e6)

# -----------------------------------------------------------------
# Extracting ancestries
# -----------------------------------------------------------------

new_gwas_info <- gwas_study_info %>%
	left_join(opengwas_dat[, c("id", "population")])

new_ewas_info <- ewas_study_info %>%
	left_join(ewascatalog_dat[, c("Ancestry", "StudyID")])

new_ewas_info %>%
	dplyr::filter(Trait %in% new_gwas_info$ec_trait) %>%
	pull(Ancestry)

new_ewas_info %>%
	dplyr::filter(Trait %in% new_gwas_info$ec_trait)

new_ewas_info[new_ewas_info$PMID == 34887417, "Ancestry"] <- "European, African"
new_ewas_info[new_ewas_info$PMID == 34887389, "Ancestry"] <- "European, African, South Asian"

## PAPERS WERE CHECKED MANUALLY BEFORE RUNNING THIS!!!!
new_gwas_info$population <- ifelse(is.na(new_gwas_info$population), "European", new_gwas_info$population)