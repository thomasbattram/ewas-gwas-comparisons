# ------------------------------------------------
# Simulations to compare methods for assessing overlap of pathways
# ------------------------------------------------

# There are 3 (kind of 2.5) methods look at whether there is more 
# overlap between pathways identified by EWAS and GWAS than expected
# by chance. This script aims to assess which is the best! 

pkgs <- c("tidyverse", "FDb.InfiniumMethylation.hg19",
		  "IlluminaHumanMethylation450kanno.ilmn12.hg19")

lapply(pkgs, require, character.only = T)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

devtools::load_all("~/repos/usefunc/")

source("R/mapping_functions.R")
source("R/simulation_functions.R")

args <- commandArgs(trailingOnly = TRUE)
split <- as.numeric(args[1])

split1 <- (split - 1) * 10 + 1
split2 <- split * 10

message("split = ", split1, " to ", split2)

# split <- as.numeric(unlist(strsplit(split, ",")))
# message("sim function = ", sim_n)
# message("split = ", split)

# gene ontology terms - table of terms + genes linked to the terms
# ensembl_gene_id pathway_id description
go_terms <- read_tsv("data/derived/gene_ontology_terms.txt") %>% 
	dplyr::filter(!is.na(pathway_id))
# same as above but for kegg
kegg_terms <- read_tsv("data/derived/kegg_terms.txt") 

# gene names + positions
# ensembl_gene_id hgnc_symbol chromosome_name start_position end_position
all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)
all_genes$size <- all_genes$end_position - all_genes$start_position
summary(all_genes$size)

# unique pathways
u_go <- unique(go_terms$pathway_id)
u_kegg <- unique(kegg_terms$pathway_id)

# ------------------------------------------------
# Setup simulations
# ------------------------------------------------

# First simulation, want to see in which situations the methods
# are able to detect true overlap above chance (i.e. when percent_ewas_causal = 0)

params <- expand.grid(
	n_genes = c(500, 1000, 2500, 5000, 10000, 25000),
	gwas_power = c(0.1),
	ewas_power = c(0.05, 0.1, 0.2, 0.5, 1), 
	percent_ewas_causal = c(0, 1),
	pathway_db = c("kegg", "go"),
	sim = c(1:1000),
	gene_overlap = NA, 
	pathway_overlap = NA, 
	or_g = NA,
	or_p = NA, 
	rho_p= NA
	)


extract_pathways <- function(gwas_genes, ewas_genes, databases)
{
	# function to extact pathways from each of the databases
	# supplied. Database dfs should be name DATABASE_terms

	out_paths <- lapply(seq_along(databases), function(x) {
		pathway_dat <- get(paste0(databases[x], "_terms"))
		gwasp <- pathway_dat %>%
			dplyr::filter(ensembl_gene_id %in% gwas_genes) %>%
			pull(pathway_id) %>%
			unique
		ewasp <- pathway_dat %>%
			dplyr::filter(ensembl_gene_id %in% ewas_genes) %>%
			pull(pathway_id) %>%
			unique
		return(list(gwas = gwasp, ewas = ewasp))
	})
	out_paths <- flatten(out_paths)
	return(out_paths)
}

extract_genes_in_pathways <- function(pathway_df) 
{
	# takes all genes from each pathway and outputs this as a list    

    unique_pathways <- unique(pathway_df$pathway_id)
    len_paths <- length(unique_pathways)
    all_pathways <- lapply(seq_along(unique_pathways), function(x) {
        message("Pathway ", x, " of ", len_paths)
    	pathway_of_interest <- unique_pathways[x]
    	pathway_genes <- pathway_df %>%
    		dplyr::filter(pathway_id %in% pathway_of_interest) %>%
    		pull(ensembl_gene_id)
    	return(pathway_genes)
    })
    names(all_pathways) <- unique_pathways
    return(all_pathways)
}

perform_enrichment <- function(identified_genes, database, background_genes) 
{
    # uses fishers exact test to perform geneset enrichment analyses 
    # for each pathway and output pathway_id and "enrichment scores"

	pathway_gene_data <- get(paste0(database, "_genes"))
	len_paths <- length(pathway_gene_data)
    enriched_paths <- map_dfr(seq_along(pathway_gene_data), function(x) {
        # message("Pathway ", x, " of ", len_paths)
        pathway_of_interest <- names(pathway_gene_data)[x]      
        pathway_genes <- pathway_gene_data[[pathway_of_interest]]

        res <- overlap_test(group1 = identified_genes, 
                            group2 = pathway_genes, 
                            all_variables = background_genes)
        out <- data.frame(pathway_id = pathway_of_interest, enrich_or = res$estimate, 
        				  enrich_p = res$p) %>%
        		mutate(adj_enrich_p = p.adjust(enrich_p, method = "BH"))
        return(out)
    }) 
}

overlap_test <- function(group1, group2, all_variables) 
{
    # asks if there is more of an overlap between
    # group one and group two than expected by chance if 
    # those groups were randomly sampled from all the variables

    q <- sum(group1 %in% group2) # overlap between ewas and gwas pathways
    m <- length(group2) - q # number of pathways identified by GWAS but not EWAS
    k <- length(group1) - q # number of pathways identified by EWAS but not by GWAS
    n <- length(all_variables) - q - m - k # all pathways minus pathways identified by EWAS and GWAS
    tab <- matrix(c(q, m, k, n), 2, 2)
    out <- fisher.test(tab, alternative = "greater")
    return(out)
}

sim_func <- function(n_genes, ewas_power, gwas_power, percent_ewas_causal, database, sim) 
{
	# simulation function

	set.seed(sim)
	pathway_dat <- get(paste0(database, "_terms"))
	# sample genes
	causal_trait_genes <- sample(all_unique_genes, n_genes)
	gwas_genes <- sample(causal_trait_genes, n_genes * gwas_power)
	n_ewas_genes <- n_genes * ewas_power
	n_ewas_causal_genes <- n_ewas_genes * percent_ewas_causal
	ewas_genes <- c(sample(causal_trait_genes, n_ewas_causal_genes), 
					sample(all_unique_genes, n_ewas_genes - n_ewas_causal_genes))
	# perform gene overlap test
	message("Testing if genes overlap")
	overlap_g <- overlap_test(ewas_genes, gwas_genes, all_unique_genes)
	# extract pathways and do overlap test on them all
	all_pathways <- extract_pathways(gwas_genes, ewas_genes, database)
	gwasp <- all_pathways$gwas
	ewasp <- all_pathways$ewas
	message("Testing if pathways overlap")
	overlap_p <- overlap_test(ewasp, gwasp, unique(pathway_dat$pathway_id))
	# perform enrichment analyses
	message("Performing enrichment analyses")
	gwas_enrich <- perform_enrichment(gwas_genes, database, all_unique_genes)
	ewas_enrich <- perform_enrichment(ewas_genes, database, all_unique_genes)
	enrich_cor <- cor(log(gwas_enrich$enrich_or), log(ewas_enrich$enrich_or),
					  use = "pairwise.complete.obs",
					  method = "spearman")

	out <- list(gene_overlap = sum(ewas_genes %in% gwas_genes), 
				gene_res = overlap_g$estimate, 
				pathway_overlap = sum(ewasp %in% gwasp), 
				or_p = overlap_p$estimate, 
				rho_p = enrich_cor)
	return(out)
}

# doing this makes performing enrichment analyses much quicker
databases <- c("kegg", "go")
path <- "data/derived"
db_genes <- lapply(databases, function(db) {
	db_gene_file <- paste0(db, "_genes.rdata")
	full_path <- file.path(path, db_gene_file)
	if (file.exists(full_path)) {
		load(full_path)
		pathway_genes <- get(paste0(db, "_genes"))
	} else {
		pathway_df <- get(paste0(db, "_terms"))
		pathway_genes <- extract_genes_in_pathways(pathway_df)
		save(pathway_genes, file = full_path)
	}
	return(pathway_genes)
})
kegg_genes <- db_genes[[1]]
go_genes <- db_genes[[2]]

# some genes are present multiple times, but start and positions are the
# same so just take unique entries
all_unique_genes <- unique(all_genes$ensembl_gene_id)

test_params <- params %>%
	dplyr::filter(sim == 2, 
				  pathway_db == "kegg", 
				  n_genes == max(n_genes), 
				  ewas_power == 1, 
				  percent_ewas_causal == 1)

params <- params %>%
	dplyr::filter(between(sim, split1, split2))

# params <- test_params

start_time <- proc.time()
out <- lapply(split(params, 1:nrow(params)), function(x) {
	print(x)
	sim_res <- sim_func(n_genes = x$n_genes, 
						ewas_power = x$ewas_power, 
						gwas_power = x$gwas_power, 
						percent_ewas_causal = x$percent_ewas_causal, 
						database = x$pathway_db, 
						sim = x$sim)
	x$gene_overlap <- sim_res$gene_overlap
	x$pathway_overlap <- sim_res$pathway_overlap
	x$or_g <- sim_res$gene_res
	x$or_p <- sim_res$or_p
	x$rho_p <- sim_res$rho_p
	return(x)
})
out <- bind_rows(out)
time_taken <- proc.time() - start_time
time_taken

path <- "data/derived/sim_temp/"
nam <- paste0("method_sim", "_", split2, ".txt")
write.table(out, file = paste0(path, nam), 
			row.names = F, col.names = T, quote = F, sep = "\t")

# if (split2 == 1000) {
# 	files <- list.files(path)
# 	sim_files <- grep("method_sim", files, value = TRUE)
# 	sim_res <- 	map_dfr(sim_files, function(fil) {
# 		x <- read_tsv(paste0(path, fil))
# 		return(x)
# 	})
# 	nam <- paste0("results/sim/method_sim_res.txt")
# 	write.table(sim_res, file = nam, 
# 				row.names = F, col.names = T, quote = F, sep = "\t")	
# }
