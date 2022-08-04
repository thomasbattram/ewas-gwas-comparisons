# ---------------------------------------------------------
# Script to sort the epigenetic data
# ---------------------------------------------------------

## aim of script:
# to extract data from the EWAS Catalog and merge in any external data
# that has been found that is not within the catalog.
# then to add in the genomic position information and output it all

## pkgs
library(tidyverse) # for tidy code + data
library(usefunc) # personal package of useful functions

# setwd(SCRATCH_SPACE)

## read in the partitioned data
cpg_divisions <- read_tsv("data/derived/epigenetic/divided_genome.txt", guess_max = 1e6)
# snp_divisions <- read_tsv("data/derived/genetic/genetic_data_divided_genome.txt")

## read in the ewas catalog data --- MAKE SURE THIS IS UP TO DATE!
# ewas_cat <- read_tsv("data/ewas_catalog/EWAS_Catalog_23-11-2018.txt", guess_max = 1e6)
ec_studies <- read_tsv("data/ewas_catalog/studies.txt")
ec_results <- read_tsv("data/ewas_catalog/results.txt")

## read in the external ewas data
# get any external data if needs be
# this should be extracted and sorted manually using "extract_ewas_data_manual.R"
files <- c(insulin = "data/glucose_ewas/insulin_ewas_data.RData", 
		   glucose = "data/glucose_ewas/glucose_ewas_data.RData", 
		   egfr = "data/ckd_ewas/egfr_ewas_data.RData", 
		   urate = "data/ckd_ewas/urate_ewas_data.RData") # Add R data files 
external_dat <- lapply(files, function(x) {
	if (!file.exists(x)) stop(paste("File", x, "doesn't exist"))
	out <- new_load(x)
	out$studies$PMID <- as.character(out$studies$PMID)
	return(out)
})

sort_epigen_dat <- function(epi_studies, epi_results, trait, cols_to_keep) 
{
	### function to bind ewas data with genomic annotations
	epi_studies %>%
		dplyr::filter(Trait == trait) %>%
		left_join(epi_results) %>%
		dplyr::filter(P < 1e-4) %>%
		dplyr::select(one_of(cols_to_keep)) %>%
		left_join(cpg_divisions)
}


# ---------------------------------------------------------
# sort the data
# ---------------------------------------------------------

## get all studies data from the external dat
## merge all studies and results data
if (exists("external_dat")) {
	names(external_dat) <- names(files)
	external_studies <- map_dfr(external_dat, "studies")
	external_results <- map_dfr(external_dat, "results")
	studies <- bind_rows(ec_studies, external_studies)
	results <- bind_rows(ec_results, external_results)
} else {
	studies <- ec_studies
	results <- ec_results
}

# first sort the ewas catalog data
large_n_dat <- studies %>%
	dplyr::filter(N > 4500)

large_dmp_dat <- results %>%
	dplyr::filter(StudyID %in% large_n_dat$StudyID) %>%
	dplyr::filter(P < 1e-7) %>%
	group_by(StudyID) %>%
	summarise(n_dmps = n()) %>%
	dplyr::filter(n_dmps > 10)

### Some manual stuff that should really be removed and fixed in the actual EWAS Catalog data...
large_n_dat <- large_n_dat %>%
	dplyr::filter(StudyID %in% large_dmp_dat$StudyID) %>%
	mutate(Trait = tolower(gsub("[[:space:]]", "_", Trait))) %>%
	mutate(Exposure = tolower(gsub("[[:space:]]", "_", Exposure))) %>%
	mutate(Outcome = tolower(gsub("[[:space:]]", "_", Outcome)))

large_n_dat %>%
	dplyr::select(Author, Trait, Exposure, Outcome, Analysis) %>%
	dplyr::filter(Trait == "smoking") %>%
	distinct() %>%
	as.data.frame
# for smoking the trait was labelled as smoking, but exposure was different! 
# also for smoking it tells you the more specific trait in the analysis section
# e.g. former vs never 
large_n_dat <- large_n_dat %>%
	mutate(Trait = ifelse(Author == "Joehanes R", Exposure, Trait))

## selecting just one model from studies
trait <- unique(large_n_dat$Trait)[7]
large_studies <- map_dfr(unique(large_n_dat$Trait), function(trait) {
	print(trait)
	dat <- large_n_dat %>%
		dplyr::filter(Trait == trait)

	out_dat <- dat %>%
		dplyr::filter(N == max(N)) %>%
		dplyr::select(PMID, Author, Trait, Exposure, Outcome, Analysis, N, StudyID)

	return(out_dat)
})

## check analyses aren't bs -- i.e. if a study adjusted for batch in one
## model but not in another then need to choose batch adjusted one!
large_studies$Analysis

## Manual changes:

## Removing pmid=31536415 as that just used data already used in 2 other EWAS
large_studies <- large_studies %>%
	dplyr::filter(PMID != 31536415)

## Removing models of certain results that say they're not adjusted for cell counts
studs_to_rm <- c("31230546_Kazmi-N_hypertensive_disorders_of_pregnancy_meta_analysis_unadjusted_model", 
				 "29016858_Sharp-GC_maternal_body_mass_index_meta-analysis", 
				 "29016858_Sharp-GC_maternal_overweight_obesity_meta-analysis")

## Removing smaller bmi EWAS
studs_to_rm <- c(studs_to_rm, "29278407_Geurts-YM_bmi_discovery")

## Removing smaller alcohol consumption EWAS
studs_to_rm <- c(studs_to_rm, "31789449_Dugue-P-A_alcohol_consumption_replication")

## Removing homocysteine grs and mutation EWAS - don't want to compare an EWAS of a GWAS generated trait to the GWAS of that trait
studs_to_rm <- c(studs_to_rm, "29084233_Mandaviya-PR_homocysteine_grs_meta-analysis", 
							  "29084233_Mandaviya-PR_homocysteine__mthfr_677c_t__meta-analysis")

## Removing glucose + insulin EWAS adjusted for BMI
studs_to_rm <- c(studs_to_rm, "31197173_Liu-J_fasting_insulin_base_plus_bmi", 
							  "31197173_Liu-J_fasting_glucose_base_plus_bmi")

## Removing old CKD EWAS
studs_to_rm <- c(studs_to_rm, "29097680_Chu-A._chronic_kidney_disease_meta-analysis_of_prevalent_ckd_ewas", 
							  "29097680_Chu-A._chronic_kidney_disease_meta-analysis_of_estimated_glomerular_filtration_rate_ewas")

large_studies <- large_studies %>%
	dplyr::filter(!StudyID %in% studs_to_rm)

## NOTES:
# multiple smoking in preg - pmid = 27040690, pmid = 31536415
# pmid = 31536415 - smoking is ever v never? If so than can keep it in
# keep row 16 (pmid = 31230546, analysis=adjusted model) and delete row 15
# check glucose + insulin ewas to see whether we previously used adjusted or 
# 		unadjusted data -- not adjusted for BMI!
# keep row 28 and row 30 (maternal BMI and maternal obesity adjusted models) and remove row 27 + 29
# remove second bmi EWAS (pmid = 29278407)
# remove prevalent CKD EWAS

## changing trait names with a slash in them so can write out the files! 
large_studies$Trait <- gsub("maternal_overweight\\/obesity", "maternal_overweight_or_obesity", large_studies$Trait) 
large_studies$Trait <- gsub("fev1_\\/_fvc", "fev1_fvc_ratio", large_studies$Trait) 

## extract results with large N 
large_results <- results %>%
	dplyr::filter(StudyID %in% large_studies$StudyID)

## sort out the ewas data
cols_to_keep <- c("Author", "Trait", "Exposure", "Outcome", "CpG", "Beta", "SE", "P", "N")
lapply(large_studies$Trait, function(trait) {
	print(trait)
	out_nam <- paste0("data/derived/epigenetic/sorted_", trait, ".txt")
	if (!file.exists(out_nam)) {
		print(paste("Making epigenetic file for", trait))
		out_epi_dat <- sort_epigen_dat(epi_studies = large_studies, 
									   epi_results = large_results, 
									   trait = trait, 
									   cols_to_keep = cols_to_keep)
		write.table(out_epi_dat, file = out_nam, 
					col.names = T, row.names = F, quote = F, sep = "\t")
	}
	return(NULL)
})

## write out the traits!
write.table(large_studies$Trait, file = "data/traits.txt",
			col.names = F, row.names = F, quote = F, sep = "\n")

## write out studies info! 
write.table(large_studies, file = "data/derived/epigenetic/ewas_studies_info.tsv", 
			col.names = T, row.names = F, quote = F, sep = "\t")

