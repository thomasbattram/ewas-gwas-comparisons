# -------------------------------------------------------
# Simulating EWAS and GWAS pathway overlap
# -------------------------------------------------------

pkgs <- c("tidyverse", "FDb.InfiniumMethylation.hg19",
		  "IlluminaHumanMethylation450kanno.ilmn12.hg19")

lapply(pkgs, require, character.only = T)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

devtools::load_all("~/repos/usefunc/")

source("R/mapping_functions.R")
source("R/simulation_functions.R")

args <- commandArgs(trailingOnly = TRUE)
sim_n <- as.numeric(args[1]) # simulation number (can be 1, 2, 3 or 4)
split <- as.numeric(args[2])
split1 <- (split - 1) * 10 + 1
split2 <- split * 10
message("split = ", split1, " to ", split2)

egdb_reactome <- read_tsv("data/derived/epigraphdb_proteins_and_pathways.tsv")
egdb_ppi <- read_tsv("data/derived/epigraphdb_proteins-protein_interactions.tsv")


# gene names + positions
# ensembl_gene_id hgnc_symbol chromosome_name start_position end_position
# all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)
# all_genes$size <- all_genes$end_position - all_genes$start_position
# summary(all_genes$size)

# -------------------------------------------------------
# Setup datasets
# -------------------------------------------------------

react_pathway_names <- egdb_reactome %>%
	pull(pathway_reactome_id) %>%
	gsub("\\[", "", .) %>%
	gsub("\\]", "", .) %>%
	gsub("\\'", "", .) %>%
	strsplit(", ") %>%
	unlist() %>%
	unique

x=react_pathway_names[1]
reactome_res <- lapply(react_pathway_names, function(x) {
	proteins_in_pathway <- egdb_reactome %>%
		dplyr::filter(grepl(x, pathway_reactome_id)) %>%
		pull(uniprot_id) %>%
		unique
})
names(reactome_res) <- react_pathway_names

reactome_terms <- map_dfr(react_pathway_names, function(x) {
	tibble(pathway_id = x, 
		   protein = reactome_res[[x]])
})


uniq_proteins <- unique(egdb_ppi$protein)
p=uniq_proteins[1]
ppi_res <- lapply(uniq_proteins, function(p) {
	associated_proteins <- egdb_ppi %>%
		dplyr::filter(protein == p) %>%
		pull(assoc_protein) %>%
		unique
	associated_proteins[!is.na(associated_proteins)]
})
names(ppi_res) <- uniq_proteins

# system.time(names(ppi_res)[grep("P32455", ppi_res)])
# system.time(egdb_ppi %>% dplyr::filter(assoc_protein=="P32455") %>% pull(protein))
# system.time(names(ppi_res)[map_lgl(ppi_res, function(x) "P32455" %in% x)])

# much faster to filter from data.frame than to search list

ppi_terms <- egdb_ppi %>%
	dplyr::select(pathway_id = protein, protein = assoc_protein) %>%
	dplyr::filter(!is.na(protein))

# -------------------------------------------------------
# Setup parameters
# -------------------------------------------------------

# # pathway down approach
# params1 <- expand.grid(
# 	prop_causal_pathways = c(0.01),
# 	prop_consequent_pathways = c(0.005, 0.01, 0.05, 0.1),
# 	gwas_power = c(0.1),
# 	ewas_power = c(0.05, 0.1, 0.2, 0.5, 1),
# 	percent_ewas_causal = c(0, 0.05, 0.1, 0.2, 0.5, 1),
# 	pathway_type = c("reactome", "ppi"),
# 	sim = c(1:1000),
# 	n_gwas_genes = NA,
# 	n_ewas_genes = NA, 
# 	gene_overlap = NA, 
# 	pathway_overlap = NA, 
# 	or_g = NA,
# 	or_p = NA, 
# 	rho_p= NA
# )

# gene up approach
params2 <- expand.grid(
	n_genes = c(500, 1000, 2500, 5000, 10000, 25000),
	gwas_power = c(0.1),
	ewas_power = c(0.05, 0.1, 0.2, 0.5, 1),
	percent_ewas_causal = c(0, 0.05, 0.1, 0.2, 0.5, 1),
	pathway_type = c("reactome", "ppi"),
	sim = c(1:1000),
	gene_overlap = NA, 
	pathway_overlap = NA, 
	or_g = NA,
	or_p = NA, 
	rho_p= NA
)

params <- get(paste0("params", sim_n))
params_to_rm <- grep("params[0-9]", ls(), value = TRUE)
rm(list = params_to_rm)

params <- params %>%
	dplyr::filter(between(sim, split1, split2))

test_params <- params %>%
	dplyr::filter(sim == 1)

# params <- test_params

# have messages set up!
message_at <- calc_increments(nrow(params))

print("starting simulations")

function_name <- paste0("apply_sim", sim_n)
sim_func <- match.fun(function_name)
start_time <- proc.time()
out <- sim_func(params)
out <- bind_rows(out)
time_taken <- proc.time() - start_time

path <- "data/derived/sim_temp/"
nam <- paste0("sim", sim_n, "_", split, ".txt")
write.table(out, file = paste0(path, nam), 
			row.names = F, col.names = T, quote = F, sep = "\t")


