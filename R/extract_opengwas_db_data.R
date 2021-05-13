# ----------------------------------------------------
# extracting top GWAS hits from opengwas project
# ----------------------------------------------------

## script uses opengwas project database to extract all gwas data
## under certain specs. This if for the analysis to check if there is 
## overlap between all gwas and the extracted ewas in the overlapping_many_gwas_signals.R
## script.

## NB. This needs to be run in BC4! 

pkgs <- c("tidyverse", "gwasvcf", "ieugwasr")
lapply(pkgs, require, character.only = T)

# devtools::load_all("~/repos/usefunc")

dir <- "~/projects/epi_gen_comp"
setwd(dir)

gwasdb_path <- "/mnt/storage/private/mrcieu/data/IGD/data/public"

# ---------------------------------------------
# get ids that are needed
# ---------------------------------------------
gwasdb_info <- gwasinfo()

selected_gwas <- gwasdb_info %>%
	arrange(desc(sample_size)) %>%
	dplyr::filter(sample_size > 5000) %>%
	dplyr::filter(population == "European") %>%
	dplyr::filter(is.na(ncase) | (ncase > 0.1*sample_size & ncase > 100)) %>%
	dplyr::filter(is.na(ncontrol) | (ncontrol > 0.1*sample_size & ncontrol > 100)) %>%
	dplyr::select(trait, id, ncase, ncontrol, sample_size, pmid)

batches <- c("bbj-a", "ebi-a", "eqtl-a", "finn-a", "ieu-a", "ieu-b",
			 "met-a", "met-b", "met-c", "prot-a", "prot-b", "prot-c", 
			 "ubm-a", "ukb-b", "ukb-d")

batches_to_rm <- "eqtl-a"



selected_gwas <- map_dfr(batches, function(b) {
	dat <- selected_gwas %>%
		dplyr::filter(grepl(b, id)) %>%
		mutate(batch = b)
	return(dat)
})

selected_gwas <- selected_gwas %>%
	dplyr::filter(!batch %in% batches_to_rm) %>%
	mutate(vcffile = paste0(id, ".vcf.gz"))

table(selected_gwas$batch) # no eqtl or protein data --> Good!

selected_gwas %>%
	group_by(batch) %>%
	summarise(min_samp = min(sample_size), min_case = min(ncase, na.rm = T), 
			  min_ncontrol = min(ncontrol, na.rm = T))


# testing which IDs are there! 
vcffile_paths <- file.path(gwasdb_path, selected_gwas$id, selected_gwas$vcffile)

x <- map_lgl(vcffile_paths, file.exists)

bad_ids <- selected_gwas[!x, "id", drop=T]

message("These ID doesn't have a corresponding vcf file: ", paste(bad_ids, collapse = " "))

writeLines(bad_ids, con = "data/bad_ids.txt")

selected_gwas <- selected_gwas[x, ]

# minimum sample size for all is above 

set.seed(2)
test_out <- selected_gwas[sample(1:nrow(selected_gwas), 300), ]

write.table(selected_gwas, file = "data/ieugwasdb_studies.txt", 
			row.names = F, col.names = T, quote = F, sep = "\t")

# ---------------------------------------------
# extract data
# ---------------------------------------------

set_bcftools(path = genetics.binaRies::get_bcftools_binary())
dp <- "data"
id <- selected_gwas[1, "id", drop = T]

extract_data <- function(id, trait, gwasdb_path, data_path = dp, 
						 tmp_path = "tmp", pthreshold = 5e-8) 
{
	### Extract number of gwas hits from ukb data in 
	### ieu open gwas database vcf files
	if (!file.exists(data_path)) stop("Data path doesn't exist")
	out_nam <- paste0(id, ".tsv")
	out_path <- file.path(data_path, tmp_path, out_nam)
	if (file.exists(out_path)) return(NULL)
	message("Extracting data for ", id)
	filename <- paste0(id, ".vcf.gz")
	vcffile <- file.path(gwasdb_path, id, filename)
	vcf <- query_gwas(vcffile, pval=pthreshold)
	dat <- vcf_to_tibble(vcf)
	if (nrow(dat) == 0) {
		out_res <- data.frame(id = id, trait = trait, SNP = "none")
	} else {
		out_res <- dat %>%
			mutate(trait = trait) %>%
			dplyr::select(id, trait, SNP = rsid)
	}

	write.table(out_res, file = out_path, 
				col.names = T, row.names = F, quote = F, sep = "\t")
}
row=1
start_time <- proc.time()
lapply(1:nrow(selected_gwas), function(row) {
	id <- selected_gwas[row, "id", drop = TRUE]
	trait <- selected_gwas[row, "trait", drop = TRUE]
	extract_data(id, trait, gwasdb_path)
})
time_taken <- proc.time() - start_time
time_taken # for 5 traits

((time_taken[["elapsed"]] / 5) * nrow(selected_gwas))/60/60
# 48h...

# ---------------------------------------------
# read in data and bind it together
# ---------------------------------------------

bind_data <- function(filepath = file.path(dp, "tmp")) 
{
	### Read in the extracted data and bind it together
	files <- list.files(filepath)
	map_dfr(files, function(file) {
		read_tsv(file.path(filepath, file))
	})
}

all_data <- bind_data()

write.table(all_data, file = file.path(dp, "ieugwasdb_top_hits.txt"), 
			col.names = T, row.names = F, quote = F, sep = "\t")

message("FIN")
