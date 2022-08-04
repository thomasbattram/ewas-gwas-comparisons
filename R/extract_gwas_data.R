# ----------------------------------------------------
# extracting GWAS from open gwas project
# ----------------------------------------------------

## Aim of script:
## Take the traits identified to have EWAS of N > 4500 + extract 
## the corresponding GWAS data

## pkgs
library(tidyverse) # tidy code + data
library(ieugwasr) # extracting data from open gwas database
library(gwasvcf) # extracting data from open gwas database on bc4 (if api doesn't work for ieugwasr)
library(usefunc) # personal package with useful functions in
library(Rmpfr) # for calculating p values

## read in traits identified from ewas catalog!
traits <- readLines("data/traits.txt")

# ----------------------------------------------------
# first pass search of open gwas project data
# ----------------------------------------------------

search_date_file <- "data/derived/date_of_gwas_search.RData"
if (file.exists(search_date_file)) {
	load(search_date_file)
	message("Last search date: ", gwas_search_date)
}

## assign keywords to traits to help with gwas search
# --- ALTER MANUALLY!!!
traits
labelled_traits <- tribble(
  ~keywords,            		~trait,
  "maternal_smoke",     		"maternal_smoking_in_pregnancy",
  "smoke",              		"current_versus_never_smoking",
  "smoke",              		"former_versus_never_smoking",
  "alcohol",            		"alcohol_consumption_per_day",
  "c-reactive_protein", 		"c-reactive_protein",
  "body_mass_index",    		"body_mass_index",
  "education",          		"educational_attainment",
  "school",             		"educational_attainment",
  "kidney_disease", 				"egfr",
  "egfr",										"egfr", 
  "urate",									"urate",
  "pre-eclampsia", 					"hypertensive_disorders_of_pregnancy",
  "hypertension", 					"hypertensive_disorders_of_pregnancy",
  "insulin",            		"fasting_insulin",
  "glucose",            		"fasting_glucose", 
  "blood_pressure", 				"systolic_blood_pressure",
  "blood_pressure", 				"diastolic_blood_pressure",
  "birthweight", 						"birthweight", 
  "cognition", 							"cognitive_abilities:_digit_test", 
  "cognitive", 							"cognitive_abilities:_digit_test", 
  "maternal_body_mass_index", 	"maternal_body_mass_index", 
  "maternal_bmi",						"maternal_body_mass_index", 
  "maternal_overweight",		"maternal_overweight_or_obesity", 
  "fev1", 									"fev1_fvc_ratio",
  "fvc", 										"fev1_fvc_ratio",
  "Forced expiratory volume", 	"fev1_fvc_ratio",
  "Forced vital capacity", 	"fev1_fvc_ratio",
  "fev1", 									"fev1",
  "Forced expiratory volume",	"fev1"
)

## GET RID OF DOUBLE ALCOHOL CONSUMPTION!!! AND HOMOCYSTEINE GRS!!

## get gwas studies info
gwasdb_info <- gwasinfo()
## sometimes this doesn't work because it can't connect to the API properly so it needs to be run 
## locally and then moved to the HPC
## gwasdb_info <- read_tsv("gwasdb_info-2022-05-31.tsv", guess_max = 1e6)
gwasdb_info$trait <- tolower(gsub(" ", "_", gwasdb_info$trait))

## perform a "fuzzy" search for the keywords for each trait within
## the gwas database
gwasdb_trait_info <- lapply(seq_along(traits), function(x) {
	tr <- traits[x]
	keywords <- labelled_traits %>%
		dplyr::filter(trait == tr) %>%
		pull(keywords)
	map_dfr(keywords, function(keyword) {
		gwasdb_info %>% 
			dplyr::filter(agrepl(keyword, trait)) %>%
			dplyr::select(id, author, pmid, trait, sample_size, nsnp) %>%
			dplyr::filter(sample_size > 5000) %>%
			dplyr::filter(nsnp > 1e6) %>%
			arrange(desc(sample_size))
	})
})
names(gwasdb_trait_info) <- traits

# save the date of the search
gwas_search_date <- Sys.Date()
save(gwas_search_date, file=search_date_file)

# ----------------------------------------------------
# manually choose GWAS from the "fuzzy" match search
# ----------------------------------------------------

traits

### maternal smoking in pregnancy
gwasdb_trait_info[["maternal_smoking_in_pregnancy"]]
message("No matches. Search date: ", gwas_search_date)

# there is maternal smoking AROUND birth, but not during pregnancy...
gwasdb_trait_info[["maternal_smoking_in_pregnancy"]] <- NULL

### current versus never smoking --- can't find right data (pmid:30643251)
gwasdb_trait_info[["current_versus_never_smoking"]]$trait
gwasdb_info %>%
	dplyr::filter(pmid == 30643251) %>%
	pull(trait)

## manually read in data!
# smok_init_file <- "data/gwas_sum_stats/liu_smoking_initiation_2019_sumstats.txt.gz"
smok_init_file <- "data/gwas_sum_stats/liu_smoking_initiation_2019_sumstats.txt"
smok_init <- data.table::fread(smok_init_file)
## extract study data info
gwasdb_trait_info[["current_versus_never_smoking"]] %>%
	mutate(id = NA, 
		   author = "Liu M", 
		   pmid = 30643251, 
		   trait = "ever_smoked_regularly", 
		   sample_size = getmode(smok_init$N), 
		   nsnp = nrow(smok_init)) %>%
		   distinct() -> gwasdb_trait_info[["current_versus_never_smoking"]]

## reduce data size
# important to make object name same as ewas catalog trait name
# AND to change column names to match those from the opengwas project		   
current_versus_never_smoking <- smok_init %>%
	dplyr::filter(PVALUE < 1e-5) %>%
	dplyr::select(rsid = RSID, beta = BETA, se = SE, p = PVALUE)
rm(smok_init)

### former versus never smoking
df <- gwasdb_trait_info[["former_versus_never_smoking"]]
df$trait
best_matches <- c("past_tobacco_smoking", "smoking_status:_previous")

## manually check IDs and decide which to use
df %>%
	dplyr::filter(trait %in% best_matches)
id_to_use <- "ukb-b-2134"

df_out <- df %>%
	dplyr::filter(id == id_to_use)
gwasdb_trait_info[["former_versus_never_smoking"]] <- df_out

### alcohol consumption per day -- can't find right data (pmid:30643251)
df <- gwasdb_trait_info[["alcohol_consumption_per_day"]]
df$trait
df %>%
	dplyr::filter(trait %in% "alcoholic_drinks_per_week")
id_to_use <- "ieu-b-73"

df_out <- df %>%
	dplyr::filter(id == id_to_use)
gwasdb_trait_info[["alcohol_consumption_per_day"]] <- df_out

### c-reactive protein
df <- gwasdb_trait_info[["c-reactive_protein"]]
df$trait
# all good matches! 
## manually check IDs and decide which to use
df
id_to_use <- "ieu-b-35"

df_out <- df %>%
	dplyr::filter(id == id_to_use)
gwasdb_trait_info[["c-reactive_protein"]] <- df_out


### body mass index
df <- gwasdb_trait_info[["body_mass_index"]]
df$trait
# all good matches! 
## manually check IDs and decide which to use
df
id_to_use <- "ieu-b-40"

df_out <- df %>%
	dplyr::filter(id == id_to_use)
gwasdb_trait_info[["body_mass_index"]] <- df_out

### educational attainment
df <- gwasdb_trait_info[["educational_attainment"]]
df$trait
df$trait[!grepl("medication", df$trait)]
best_matches <- c("age_completed_full_time_education", "year_ended_full_time_education", 
				  "years_of_schooling")
## manually check IDs and decide which to use
df %>%
	dplyr::filter(trait %in% best_matches)
id_to_use <- "ieu-a-1239"
new_author <- "Lee JJ" # adding this as author in opengwas proj wasn't complete

df_out <- df %>%
	dplyr::filter(id == id_to_use) %>%
	mutate(author = new_author)
gwasdb_trait_info[["educational_attainment"]] <- df_out

### egfr
df <- gwasdb_trait_info[["egfr"]]
df$trait
best_matches <- c("chronic_kidney_disease", "estimated_glomerular_filtration_rate_(egfr)", 
									"serum_creatinine_(egfrcrea)", "serum_cystatin_c_(egfrcys)")

# egfr_gwas_status <- download.file(url = "https://ckdgen.imbi.uni-freiburg.de/files/Stanzick2021/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz", 
# 								  destfile = "data/gwas_sum_stats/stanzick2021_egfr.annot.gc.gz",
# 								  method = "wget")

egfr_file <- "data/gwas_sum_stats/stanzick2021_egfr.annot.gc.gz"
egfr <- data.table::fread(egfr_file)
## extract study data info
gwasdb_trait_info[["egfr"]] %>%
	mutate(id = NA, 
		   author = "Stanzick KJ", 
		   pmid = 34272381, 
		   trait = "egfr", 
		   sample_size = getmode(egfr$n), 
		   nsnp = nrow(egfr)) %>%
		   distinct() -> gwasdb_trait_info[["egfr"]]

egfr <- egfr %>%
	dplyr::filter(P.value < 1e-5) %>%
	dplyr::select(rsid = RSID, beta = Effect, se = StdErr, p = P.value)

### urate (manual one)
# urate_gwas_status <- download.file(url = "https://ckdgen.imbi.uni-freiburg.de/files/Tin2019/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt.gz", 
# 								  destfile = "data/gwas_sum_stats/tin2019_urate.txt.gz",
# 								  method = "wget")
urate_file <- "data/gwas_sum_stats/tin2019_urate.txt.gz"
urate <- data.table::fread(urate_file, fill = TRUE)
urate$n_total_sum <- ifelse(is.na(urate$n_total_sum), urate[["P-value"]], urate$n_total_sum)
urate$n_total_sum <- as.numeric(urate$n_total_sum)
urate[["P-value"]] <- pnorm(abs(urate$Effect)/urate$StdErr, lower.tail = F) * 2
## extract study data info
gwasdb_trait_info[["urate"]] %>%
	mutate(id = NA, 
		   author = "Tin A", 
		   pmid = 31578528, 
		   trait = "urate", 
		   sample_size = getmode(urate$n_total_sum), 
		   nsnp = nrow(urate)) %>%
		   distinct() -> gwasdb_trait_info[["urate"]]

urate <- urate %>%
	dplyr::filter(`P-value` < 1e-5) %>%
	dplyr::select(rsid = RSID, beta = Effect, se = StdErr, p = `P-value`)

## Calculate p values as many are 0...
b <-  mpfr(urate$beta, precBits=1000)	
se <- mpfr(urate$se, precBits=1000)	
#
# obtain two-tailed p-value
#
p <- pnorm(abs(b)/se, lower.tail=F)*2
#
# format the high precision number with 4-digit precision
#
urate_actual_p <- formatMpfr(p, digits=4)
urate$p <- as.numeric(urate_actual_p)
rm(p)

### hypertensive disorders of pregnancy
df <- gwasdb_trait_info[["hypertensive_disorders_of_pregnancy"]]
df$trait

# NOTE: Putting null because the number of cases is only ~1800 for the matching GWAS
gwasdb_trait_info[["hypertensive_disorders_of_pregnancy"]] <- NULL

### insulin
df <- gwasdb_trait_info[["fasting_insulin"]]
df$trait
best_matches <- c("fasting_blood_insulin")
## manually check IDs and decide which to use
df %>%
	dplyr::filter(trait %in% best_matches)
id_to_use <- "ebi-a-GCST005185"

df_out <- df %>%
	dplyr::filter(id == id_to_use)
gwasdb_trait_info[["fasting_insulin"]] <- df_out

### glucose
df <- gwasdb_trait_info[["fasting_glucose"]]
df$trait
best_matches <- c("glucose", "fasting_blood_glucose")
## manually check IDs and decide which to use
df %>%
	dplyr::filter(trait %in% best_matches)
id_to_use <- "ebi-a-GCST005186"

df_out <- df %>%
	dplyr::filter(id == id_to_use)
gwasdb_trait_info[["fasting_glucose"]] <- df_out

### systolic blood pressure
df <- gwasdb_trait_info[["systolic_blood_pressure"]]
df$trait
best_matches <- c("systolic_blood_pressure")
## manually check IDs and decide which to use
df %>%
	dplyr::filter(trait %in% best_matches)
id_to_use <- "ieu-b-38"

df_out <- df %>%
	dplyr::filter(id == id_to_use)
gwasdb_trait_info[["systolic_blood_pressure"]] <- df_out

### diastolic blood pressure
df <- gwasdb_trait_info[["diastolic_blood_pressure"]]
df$trait
best_matches <- c("diastolic_blood_pressure")
## manually check IDs and decide which to use
df %>%
	dplyr::filter(trait %in% best_matches)
id_to_use <- "ieu-b-39"

df_out <- df %>%
	dplyr::filter(id == id_to_use)
gwasdb_trait_info[["diastolic_blood_pressure"]] <- df_out

### birthweight
df <- gwasdb_trait_info[["birthweight"]]
df$trait
best_matches <- c("birth_weight")
## manually check IDs and decide which to use
df %>%
	dplyr::filter(trait %in% best_matches)
id_to_use <- "ieu-a-1083"
new_author <- "Horikoshi M"

# NOTE: Could have used UKB GWAS as it has a larger N, but it has a smaller number of snps tested.
#		Plus I don't trust the results as much

df_out <- df %>%
	dplyr::filter(id == id_to_use) %>%
	distinct() %>%
	mutate(author = new_author)
gwasdb_trait_info[["birthweight"]] <- df_out

### cognitive abilities
df <- gwasdb_trait_info[["cognitive_abilities:_digit_test"]]
df$trait
best_matches <- c("cognitive_performance")
## manually check IDs and decide which to use
df %>%
	dplyr::filter(trait %in% best_matches)
id_to_use <- "ebi-a-GCST006572"

# NOTE: cognitive performance (as far as I can tell) is a composite measure of many cognitive
#		measures, including the digit test

df_out <- df %>%
	dplyr::filter(id == id_to_use)
gwasdb_trait_info[["cognitive_abilities:_digit_test"]] <- df_out

### maternal body mass index
df <- gwasdb_trait_info[["maternal_body_mass_index"]]
df$trait

gwasdb_trait_info[["maternal_body_mass_index"]] <- NULL

### maternal overweight or obesity
df <- gwasdb_trait_info[["maternal_overweight_or_obesity"]]
df$trait

gwasdb_trait_info[["maternal_overweight_or_obesity"]] <- NULL

### fev1
df <- gwasdb_trait_info[["fev1"]]
df$trait
best_matches <- c("forced_expiratory_volume_in_1-second_(fev1)")
## manually check IDs and decide which to use
df %>%
	dplyr::filter(trait %in% best_matches)
id_to_use <- "ukb-b-19657"

df_out <- df %>%
	dplyr::filter(id == id_to_use) %>%
	distinct()
gwasdb_trait_info[["fev1"]] <- df_out

### fev1/fvc
df <- gwasdb_trait_info[["fev1_fvc_ratio"]]
df$trait

## NADA

gwasdb_trait_info[["fev1_fvc_ratio"]] <- NULL


# ----------------------------------------------------
# Extract the results files! 
# ----------------------------------------------------

## read in divided genome data 
snp_divisions <- data.table::fread("data/derived/genetic/genetic_data_divided_genome.txt")

## extract traits of interest
all_info <- bind_rows(gwasdb_trait_info, .id="ec_trait")

archive_old_data <- function(file_nam, main_dir)
{
	### archiving old gwas data -- needs to be done as there is 
	### some manual sorting of the data so won't be able to 
	### re-create results if data isn't archived!!
	# make archive directory
	the_date <- Sys.Date()
	dir_nam <- paste0(main_dir, "/archived_data_", the_date)
	if (!file.exists(dir_nam)) system(paste("mkdir", dir_nam))
	# move file into archived data
	system(paste("mv", file_nam, dir_nam))
	message("Archived ", file_nam)
}

x <- 1
out_dir <- "data/derived/genetic"

message("Testing we can connect to the IEU Open GWAS Project API")
api_test <- try(api_status())
lapply(1:nrow(all_info), function(x) {
	print(x)
	df <- all_info[x, ]
	id <- df$id
	print(id)
	## extract gwas hits
	if (is.na(id)) {
		# for when gwas have been manually read in
		gwas_res <- get(df$ec_trait)
	} else {
		if (inherits(api_test, "try-error")) {
			message("To perform this step you must have access to BC4")
			gwasdb_path <- "/mnt/storage/private/mrcieu/data/IGD/data/public"
			filename <- paste0(id, ".vcf.gz")
        	vcffile <- file.path(gwasdb_path, id, filename)
			gwas_res <- gwasvcf::query_gwas(vcffile, pval = 1e-5) %>% 
				gwasvcf::vcf_to_tibble() %>%
				mutate(p = 10^(-LP)) %>%
				dplyr::select(rsid, beta = ES, se = SE, p)
		} else {
			gwas_res <- tophits(id, 
								pval = 1e-5, 
								clump = 0, 
								access_token = check_access_token()
								)
		}
	}	
	message("gwas results for ", df$trait, " extracted and it is ", nrow(gwas_res), " rows")
	## sort out data
	out <- gwas_res %>%
		dplyr::select(SNP = rsid, b = beta, se, p) %>%
		left_join(snp_divisions, by = c("SNP" = "ID")) %>%
		dplyr::select(SNP, CHR, bp, b, se, p, group, group_total, positions) %>%
		dplyr::filter(!is.na(CHR))

	## write out data
	out_nam <- file.path(out_dir, paste0("sorted_", df$ec_trait, ".txt"))
	if (file.exists(out_nam)) {
		archive_old_data(out_nam, out_dir)
	}
	write.table(out, file = out_nam, 
				col.names = T, row.names = F, quote = F, sep = "\t")
})

# ----------------------------------------------------
# Write out studies info for paper
# ----------------------------------------------------

## Format so author names are correct! 
# remove commas
all_info$author <- gsub(",", "", all_info$author)
# change Ben's name to right formatting
all_info$author <- gsub("Ben Elsworth", "Elsworth B", all_info$author)

write.table(all_info, file = file.path(out_dir, "gwas_studies_info.tsv"), 
			col.names = T, row.names = F, quote = F, sep = "\t")

## write out the list of traits that you could get GWAS data
write.table(all_info$ec_trait, file = "data/traits.txt", 
			col.names = F, row.names = F, quote = F, sep = "\n")
