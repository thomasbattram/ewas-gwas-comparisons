# -------------------------------------------------------
# which gwas signal is most similar to the big gwas?
# -------------------------------------------------------

# ieu-a-1007
# ebi-a-GCST006327

pkgs <- c("tidyverse", "cluster", "ggdendro", "gplots", "usefunc")
lapply(pkgs, require, character.only = T)

args <- commandArgs(trailingOnly = TRUE)
path_type <- args[1]
split <- as.numeric(args[2])
# path_type <- "go"
# split <- 1
split1 <- (split - 1) * 10
split2 <- split * 10
message("The pathway is ", path_type)
message("split = ", split1, " to ", split2)

stopifnot(path_type %in% c("kegg", "go"))

dir <- "~/main_project/epi_gen_comp"
setwd(dir)

source("R/mapping_functions.R")
source("R/simulation_functions.R")
# check_all_dirs(trait)
# test_nam <- mrbase_snp_data.txt
# nam <- mrbase_gwas_top_hits.txt
# mrbase_snps <- read_tsv("data/derived/mrbase_gwas_top_hits.txt")
# mrbase_outcomes <- read_tsv("data/derived/mrbase_extracted_outcomes.txt")
opengwas_snps <- read_tsv("data/derived/ieugwasdb_top_hits.txt")
traits <- readLines("data/traits.txt")

# split up res file so can pararellize submissions to bc3
unique_ids <- unique(opengwas_snps$id)
task_length <- length(unique_ids)

res_file_split <- calc_increments(task_length, percent_increments = 10)

if (task_length > max(res_file_split)) {
	res_file_split['100'] <- task_length
}

res_file_split['0'] <- 0

res_file_split1 <- res_file_split[as.character(split1)] + 1
res_file_split2 <- res_file_split[as.character(split2)] 

ids_to_use <- unique_ids[res_file_split1:res_file_split2]

#ids_to_use <- 'ebi-a-GCST006327'

opengwas_snps <- opengwas_snps %>%
	dplyr::filter(id %in% ids_to_use)

opengwas_studies <- read_tsv("data/derived/ieugwasdb_studies.txt")
snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_genes.txt")
cpg_genes <- read_tsv("data/derived/epigenetic/cpgs_mapped_to_genes.txt")
pathway_dat_nam <- paste0("data/derived/", path_type, "_terms.txt")
pathway_dat <- read_tsv(pathway_dat_nam)
all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)

# map snps to genes
res <- opengwas_snps %>%
	dplyr::filter(SNP != "none") %>%
	left_join(snp_genes, by = c("SNP" = "name")) %>%
	dplyr::select(-hgnc_symbol)

ids_to_use <- unique(res$id)

print(object.size(res), units = "auto", standard = "SI")
# will probs have to split by trait before joining...
# ---------------------------------------------------------
# do some enrichment tests!!!!
# ---------------------------------------------------------

# want a matrix with each column being a trait and rows being enrichment scores
# then can just use cor function
background_snp_genes <- unique(snp_genes$ensembl_gene_id)
rm(snp_genes)

# extract only the go terms needed --> for speed! 
path <- "data/derived"
db_gene_file <- paste0(path_type, "_genes.rdata")
full_path <- file.path(path, db_gene_file)
if (file.exists(full_path)) {
    load(full_path)
    assign(paste0(path_type, "_genes"), pathway_genes)
} else {
    pathway_df <- get(paste0(db, "_terms"))
    pathway_genes <- extract_genes_in_pathways(pathway_df)
    save(pathway_genes, file = full_path)
    assign(paste0(path_type, "_genes"), pathway_genes)
}

# unique_ids <- unique_ids[1:2]
print(length(ids_to_use))
start_time <- proc.time()
i=ids_to_use[1]
g_res <- map_dfc(ids_to_use, function(i) {
	print(i)
	# pull genes of interest
	genes_of_interest <- res %>%
		dplyr::filter(id == i) %>%
		pull(ensembl_gene_id) %>%
		unique
	# generate enrichment scores for each pathway
	enrich_res <- perform_enrichment(identified_genes = genes_of_interest, 
                                  	 database = path_type, 
                                  	 background_genes = all_genes$ensembl_gene_id)
	# sanity check! 
	if (nrow(enrich_res) < length(pathway_genes)) stop("Length of results should be equal to the number of pathways")
	# output after arranging pathway ids so all should be in the same order! 
	out <- enrich_res %>%
		arrange(pathway_id) %>%
		dplyr::select(enrich_or)
	colnames(out) <- i
	return(out)
})
time_taken <- proc.time() - start_time
time_taken
ids_to_match_rows <- sort(names(pathway_genes))

nam <- paste0("gwas_data_", path_type, "_enrichment_results_", split, ".txt")
write.table(g_res, file = file.path("data/derived/enrich_res_temp", nam), 
			row.names = F, col.names = T, quote = F, sep = "\t")

# enrich_scores <- read_tsv("results/gwas_and_ewas_enrichment_scores.txt")

# x <- map_lgl(enrich_scores, function(x) sum(x > 0) == 0)
# temp <- enrich_scores[,x]

# colnames(enrich_scores) <- gsub("\\.\\.\\.[0-9]", "", colnames(enrich_scores))

# test <- dplyr::select(enrich_scores, "ieu-a-1007")
# test <- test[,,drop=T]

# test[test > 0]

# -----------------------------------------------------------------
# now do the same for the EWAS data! 
# -----------------------------------------------------------------

if (split != 1) {
	stop("Finishing analysis here as EWAS results should have already been generated")
}

start_time <- proc.time()
ewas_res <- map_dfc(traits, function(tra) {
	print(tra)
	ewas_dat <- suppressMessages(read_derived_data(tra))
	# pull genes of interest
	sig_cpgs <- ewas_dat$meth[ewas_dat$meth$P < 1e-7, "CpG", drop = T]
	genes_of_interest <- cpg_genes %>%
		dplyr::filter(name %in% sig_cpgs) %>%
		pull(ensembl_gene_id) %>%
		unique
	# generate enrichment scores for each pathway
	enrich_res <- perform_enrichment(identified_genes = genes_of_interest, 
                                  	 database = path_type, 
                                  	 background_genes = all_genes$ensembl_gene_id)
	# sanity check! 
	if (nrow(enrich_res) < length(pathway_genes)) stop("Length of results should be equal to the number of pathways")
	# output after arranging pathway ids so all should be in the same order! 
	out <- enrich_res %>%
		arrange(pathway_id) %>%
		dplyr::select(enrich_or)
	colnames(out) <- tra
	return(out)
})
time_taken <- proc.time() - start_time
time_taken

nam <- paste0("ewas_data_", path_type, "_enrichment_results_", split, ".txt")
write.table(ewas_res, file = file.path("data/derived/enrich_res_temp", nam), 
			row.names = F, col.names = T, quote = F, sep = "\t")

print("WRITTEN OUT RESULTS!!!")
