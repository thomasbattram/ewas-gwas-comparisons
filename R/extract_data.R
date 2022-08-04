# ---------------------------------------------------------
# Script to extract the genetic epigenetic data along with groups
# ---------------------------------------------------------

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = T)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

devtools::load_all("~/repos/usefunc/")

# read in the partitioned data
cpg_divisions <- read_tsv("data/derived/epigenetic/divided_genome.txt", guess_max = 1e6)
snp_divisions <- read_tsv("data/derived/genetic/genetic_data_divided_genome.txt")

# read in the ewas catalog data
ewas_cat <- read_tsv("data/ewas_catalog/EWAS_Catalog_23-11-2018.txt", guess_max = 1e6)

# read in the external ewas data
glu_ewas <- read_tsv("data/glucose_ewas/FG.model1.txt") # model 1 is without adjustment for BMI
ins_ewas <- read_tsv("data/glucose_ewas/FI.model1.txt") # model 1 is without adjustment for BMI

# genetic data files
g_dir <- "data/gwas_sum_stats/"
g_files <- paste0(g_dir, list.files(g_dir))

# traits 
traits <- c("alcohol_consumption_per_day", "former_versus_never_smoking", "current_versus_never_smoking",
		 "body_mass_index", "glucose", "insulin", "educational_attainment")
names(g_files) <- traits
variant_col_name <- c("snp", "markername", "rsid")

# ---------------------------------------------------------
# functions for sorting the genetic and epigenetic data! 
# ---------------------------------------------------------
sort_gen_dat <- function(gen_dat) {
	colnames(gen_dat) <- tolower(colnames(gen_dat))
	cols <- colnames(gen_dat)
	snp_col <- grep(paste(variant_col_name, collapse = "|"), cols, ignore.case = TRUE, value = TRUE)
	beta_col <- grep("beta|^b$|^MainEffects$", cols, ignore.case = TRUE, value = TRUE)
	se_col <- grep("^se$|^MainSE$", cols, ignore.case = TRUE, value = TRUE)
	p_col <- grep("^p$|pva.*|^MainP$", cols, ignore.case = TRUE, value = TRUE)

	new_cols <- c(snp_col, beta_col, se_col, p_col)

	if (length(new_cols) != 4) stop("Column names missing!!")

	if ("CHR" %in% cols) gen_dat <- dplyr::select(gen_dat, -CHR)
	if ("bp" %in% cols) gen_dat <- dplyr::select(gen_dat, -bp)

	out_gen_dat <- gen_dat %>%
		mutate(SNP = !! as.name(snp_col)) %>%
		mutate(b = !! as.name(beta_col)) %>%
		mutate(se = !! as.name(se_col)) %>%
		mutate(p = !! as.name(p_col)) %>%
		left_join(snp_divisions, by = c("SNP" = "ID")) %>%
		dplyr::select(SNP, CHR, bp, b, se, p, group, group_total, positions) %>%
		dplyr::filter(!is.na(CHR))

	return(out_gen_dat)
}

sort_epigen_dat <- function(epigen_dat, trait) {
	out_dat <- epigen_dat %>%
		dplyr::filter(Trait == trait) %>%
		left_join(cpg_divisions) %>%
		dplyr::filter(P < 1e-4)
}

# ----------------------------------------------
# merge the glucose and insulin ewas
# ----------------------------------------------
gi_ewas <- list(glucose = glu_ewas, insulin = ins_ewas)
var_of_interest <- c("cpg", "beta", "se", "pval", "N")

glu_dat <- map_dfr(1:2, function(x) {
	# get data
	ewas_dat <- gi_ewas[[x]] %>%
		dplyr::select(one_of(var_of_interest)) %>%
		mutate(N = max(N)) %>%
		rename(CpG = cpg, Beta = beta, SE = se, P = pval)
	ewas_nam <- names(gi_ewas)[x]
	# make data.frame that is just like the ewas catalog data.frame
	out_dat <- data.frame(matrix(ncol = length(colnames(ewas_cat)), nrow = nrow(ewas_dat)))
	colnames(out_dat) <- colnames(ewas_cat)
	cols <- colnames(ewas_dat)
	out_dat[cols] <- lapply(cols, function(i) {
		return(ewas_dat[[i]])
	})
	out_dat$Author <- "Liu J"
	out_dat$Trait <- ewas_nam
	out_dat$Exposure <- ewas_nam
	out_dat$Outcome <- "DNA methylation"
	out_dat$PMID <- 31197173
	return(out_dat)
})

merge_dat <- bind_rows(list(ewas_cat, glu_dat))

# ---------------------------------------------------------
# sort the data
# ---------------------------------------------------------

# first sort the ewas catalog data
large_n_dat <- merge_dat %>%
	dplyr::filter(N > 4500) %>%
	group_by(Trait) %>%
	summarise(count = n()) %>%
	dplyr::filter(count > 50)

large_n_dat <- merge_dat %>%
	dplyr::filter(Trait %in% large_n_dat$Trait) %>%
	dplyr::filter(N > 4500) %>%
	mutate(Trait = gsub("Alcohol consumption der day", "Alcohol consumption per day", Trait)) %>%
	mutate(Exposure = gsub("Alcohol consumption der day", "Alcohol consumption per day", Exposure)) %>%
	mutate(Trait = tolower(gsub("[[:space:]]", "_", Trait))) %>%
	mutate(Exposure = gsub("DNA mehtylation", "DNA methylation", Exposure)) %>%
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

large_dat <- map_dfr(unique(large_n_dat$Trait), function(trait) {
	print(trait)
	dat <- large_n_dat %>%
		dplyr::filter(Trait == trait)

	out_dat <- dat %>%
		dplyr::filter(N == max(N)) %>%
		dplyr::select(PMID, Author, Trait, Exposure, Outcome, CpG, Beta, SE, P, N)

	return(out_dat)
})

lapply(traits, function(tr) {
	print(tr)
	nam <- paste0("data/derived/genetic/sorted_", gsub("[[:space:]]", "_", tr), ".txt")
	if (!file.exists(nam)) {
		print(paste("Making genetic file for", tr))
		genetic_data <- read_tsv(g_files[[tr]])
		out_gen_dat <- sort_gen_dat(gen_dat = genetic_data)
		write.table(out_gen_dat, file = nam, col.names = T, row.names = F, quote = F, sep = "\t")
	}
	nam_epi <- paste0("data/derived/epigenetic/sorted_", gsub("[[:space:]]", "_", tr), ".txt")
	if (!file.exists(nam_epi)) {
		print(paste("Making epigenetic file for", tr))
		out_epi_dat <- sort_epigen_dat(epigen_dat = large_dat, trait = tr)
		write.table(out_epi_dat, file = nam_epi, col.names = T, row.names = F, quote = F, sep = "\t")
	}
	return(NULL)
})

# make table describing data
ewas_studies <- large_dat %>%
	group_by(Trait, Author, PMID) %>%
	summarise(`ewas-n` = getmode(N)) %>%
	rename(`ewas-author`= Author, `ewas-pmid` = PMID) %>%
	dplyr::filter(Trait %in% traits) 

gwas_studies <- tibble(Trait = traits, 
					   `gwas-author` = c("Liu M", "Liu M", "Liu M", "Locke AE", "Manning AK", 
					   			  		 "Manning AK", "Okbay A"),
					   `gwas-pmid` = c("30643251", "30643251", "30643251", "25673413", "22581228", "22581228",
					   				   "27225129"),  
					   `gwas-n` = c(537349, 632802, 500000, 339224, 58074, 51750, 405072))

out_study_dat <- ewas_studies %>%
	left_join(gwas_studies) %>%
	as.data.frame() %>%
	mutate(trait = gsub("_", " ", Trait)) %>%
	dplyr::select(-Trait) %>%
	dplyr::select(trait, everything())

write.table(out_study_dat, file = "report/report_data/study_data.txt", 
			col.names = T, row.names = F, quote = F, sep = "\t")
