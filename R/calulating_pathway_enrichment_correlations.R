# --------------------------------------------------
# calculating correlations between enrichment scores
# --------------------------------------------------

pkgs <- c("tidyverse", "gplots", "ggrepel", "ggExtra", "gridExtra", "usefunc")
lapply(pkgs, require, character.only = TRUE)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

args <- commandArgs(trailingOnly = TRUE)
# args <- 'go'
pathway_db <- args[1]

# --------------------------------------------------
# bind data together
# --------------------------------------------------

bind_results <- function(data_path, res_dir) 
{
	### read in and bind results from simulations
	out_nam <- file.path(res_dir, "gwas_and_ewas_enrichment_scores.txt")
	# if (file.exists(out_nam)) {
	# 	cont <- readline("File exists already. Do you want to continue? [y/n] ")
	# 	if (cont == "n") return(NULL)
	# }
	all_files <- list.files(data_path)
	message("Successful iterations: ", length(all_files))
	out <- map_dfc(all_files, function(fi) {
		suppressMessages(dat <- read_tsv(file.path(data_path, fi)))
		message("returning data")
		return(dat)
	})
	out_nam <- file.path(res_dir, "gwas_and_ewas_enrichment_scores.txt")
	message("Saving results to ", out_nam)
	write.table(out, out_nam, 
				col.names = T, row.names = F, quote = F, sep = "\t")
	return(NULL)
}

dp <- "data/derived/enrich_res_temp"
rd <- "results"

bind_results(dp, rd)


# --------------------------------------------------
# read in results and stuff
# --------------------------------------------------

# opengwas_snps <- read_tsv("data/derived/ieugwasdb_top_hits.txt")
# snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_genes.txt")

# snp_genes <- snp_genes %>%
# 	dplyr::filter(name %in% opengwas_snps$SNP) %>%
# 	dplyr::select(-chromosome, -position, -hgnc_symbol, -distance_to_gene)

# opengwas_snps <- opengwas_snps %>%
# 	dplyr::filter(SNP != "none") %>%
# 	left_join(snp_genes, by = c("SNP" = "name"))

opengwas_studies <- read_tsv("data/derived/ieugwasdb_studies.txt")
enrich_scores <- read_tsv("results/gwas_and_ewas_enrichment_scores.txt")
sum(duplicated(colnames(enrich_scores))) # 0, ayyyyyy

# remove values without any non-zeros

all_zero <- map_lgl(enrich_scores, function(x) all(x == 0))
enrich_scores <- enrich_scores[,!all_zero]

count_na <- function(dat, col_or_row = 2) 
{
	stopifnot(col_or_row %in% c(1,2))
	x <- apply(dat, col_or_row, function(x) {sum(is.na(x))})
	return(x)
}

## do correlation tests
cor_test <- function(dat) 
{
	cor_traits <- combn(colnames(dat), 2)
	message_vals <- round(seq(1, ncol(cor_traits), len = 21))
	out_res <- map_dfr(1:ncol(cor_traits), function(x) {
		# if (x %in% message_vals) {
		# 	percent_comp <- 5 * (which(message_vals == x) - 1)
		# 	message(percent_comp, "% complete")
		# }
		tra1 <- cor_traits[1,x]
		tra2 <- cor_traits[2,x]
		cor_res <- cor.test(log(dat[, tra1, drop=T]), log(dat[, tra2, drop=T]), 
							use = "pairwise.complete.obs", method = "spearman")
		if (is.na(cor_res$estimate)) return(NULL)
		out <- tibble(trait1 = tra1, 
					  trait2 = tra2, 
					  rho = cor_res$estimate, 
					  p = cor_res$p.value)
		return(out)
	})
	return(out_res)
}

# ------------------------------------------------------
# look at correlation between same GWAS traits
# ------------------------------------------------------

dup_traits <- unique(opengwas_studies[duplicated(opengwas_studies$trait), "trait", drop=T])

dup_studies <- opengwas_studies %>%
	dplyr::filter(trait %in% dup_traits) %>%
	dplyr::filter(id %in% colnames(enrich_scores))

dup_traits <- unique(dup_studies$trait)

## table wanted:
# trait id1 id2 cor cor_p
#  
#

sum(dup_studies$id %in% colnames(enrich_scores))

dt <- dup_traits[1]
dt <- 'Neuroticism'
duplicated_cor_res <- map_dfr(dup_traits, function(dt) {
	print(dt)
	study_df <- dup_studies %>%
		dplyr::filter(trait == dt)
	enrich_dat <- enrich_scores[, colnames(enrich_scores) %in% study_df$id]
	if (ncol(enrich_dat) == 0 | ncol(enrich_dat) == 1) return(NULL)
	cor_out <- cor_test(enrich_dat)
	if (nrow(cor_out) == 0) return(NULL)
	out <- cor_out %>%
		mutate(trait = dt) %>%
		dplyr::select(trait, id1 = trait1, id2 = trait2, rho, p)
	return(out)
}) 

write.table(duplicated_cor_res, file = "data/derived/duplicated_gwas_enrichment_score_correlations.txt", 
			col.names = T, row.names = F, quote = F, sep = "\t")

# ------------------------------------------------------
# generate correlations and make heatmaps
# ------------------------------------------------------
# test_enrich_scores <- enrich_scores
# remove duplicated traits! 
dup_ids_to_rm <- dup_studies %>%
	arrange(desc(sample_size)) %>%
	dplyr::filter(!duplicated(trait)) %>%
	pull(id)

enrich_scores <- enrich_scores %>%
	dplyr::select(-one_of(dup_ids_to_rm))

dim(enrich_scores)
x <- seq(10, 100, 10)
lapply(x, function(i) {
	start_time <- proc.time()
	test_enrich_scores <- enrich_scores[, 1:i]
	correlations <- cor(test_enrich_scores, use = "pairwise.complete.obs", method = "spearman")	
	time_taken <- proc.time() - start_time
	print(time_taken)
})

start_time <- proc.time()
correlations <- cor(enrich_scores, use = "pairwise.complete.obs", method = "spearman")
time_taken <- proc.time() - start_time
time_taken ### should take just over an hour to run!

# remove traits that are NA 
na_num <- count_na(correlations)
summary(na_num)
bad_traits <- names(na_num)[na_num == max(na_num)]
if (max(na_num) == 0) bad_traits <- "" 
en_dat <- enrich_scores[, !colnames(enrich_scores) %in% bad_traits]
correlations <- correlations[!rownames(correlations) %in% bad_traits, !colnames(correlations) %in% bad_traits]

save(correlations, file = "data/derived/all_enrichment_score_correlations.RData")

heatnam <- paste0("results/plots/", pathway_db, "_all_correlations_heatmap.pdf")
pdf(heatnam)
heatmap.2(as.matrix(correlations), trace = "none", scale = "none", cexRow = 0.000001, cexCol = 0.000001)
dev.off()

heatnam <- paste0("report/report_data/figure_dat/", pathway_db, "_all_correlations_heatmap.png")
png(heatnam)
heatmap.2(as.matrix(correlations), trace = "none", scale = "none", cexRow = 0.000001, cexCol = 0.000001)
dev.off()

# just ewas against everything!
ewas_traits <- readLines("data/traits.txt")
ewas_correlations <- correlations[rownames(correlations) %in% ewas_traits, !colnames(correlations) %in% ewas_traits]
heatnam <- paste0("results/plots/", pathway_db, "_ewas_correlations_heatmap.pdf")
pdf(heatnam)
heatmap.2(as.matrix(ewas_correlations), trace = "none", scale = "none", cexCol = 0.000001)
dev.off()

heatnam <- paste0("report/report_data/figure_dat/", pathway_db, "_ewas_correlations_heatmap.png")
png(heatnam)
heatmap.2(as.matrix(ewas_correlations), trace = "none", scale = "none", cexCol = 0.000001)
dev.off()

# ------------------------------------------------------
# compare correlations
# ------------------------------------------------------
start_time <- proc.time()
cor_test_res <- cor_test(en_dat)
end_time <- proc.time() - start_time # should take ~1.5h
end_time

# save all data
write.table(cor_test_res, file = paste0("results/", pathway_db, "_all_pathway_enrichment_correlations.txt"),
			col.names = T, row.names = F, quote = F, sep = "\t")


