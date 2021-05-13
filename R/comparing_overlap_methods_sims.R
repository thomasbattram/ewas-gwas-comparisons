# -------------------------------------------------------
# Testing methods to detect pathway overlap
# -------------------------------------------------------

## This script provides simulations for assessing simple methods used to
## assess overlap between genes and pathways of corresponding EWAS/GWAS

pkgs <- c("tidyverse", "FDb.InfiniumMethylation.hg19",
		  "IlluminaHumanMethylation450kanno.ilmn12.hg19")

lapply(pkgs, require, character.only = T)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

devtools::load_all("~/repos/usefunc/")

source("R/mapping_functions.R")
source("R/simulation_functions.R")

args <- commandArgs(trailingOnly = TRUE)
# args <- c("2", "1", "proteins", "go, kegg, ppi, reactome")
sim_n <- as.numeric(args[1]) # simulation number can be 1 or 2
split <- as.numeric(args[2])
genes_or_proteins <- args[3]
database_arg <- args[4]
split1 <- (split - 1) * 10 + 1
split2 <- split * 10

if (sim_n == 1) message("Pathway down simulations") ## don't use sim1
if (sim_n == 2) message("Gene up simulations")
if (!sim_n %in% 1:2) stop("sim_n has to be either 1 or 2")
message("split = ", split1, " to ", split2)
# message("sim number is going to be: ", split)

if (!genes_or_proteins %in% c("genes", "proteins")) stop("genes or proteins (argument 3) has to be either 'genes' or 'proteins'. DUH!")
message("Using the mapping of genome regions to: ", genes_or_proteins)
databases <- trimws(unlist(strsplit(database_arg, ",")))

message("Databases being used are: ", paste(databases, collapse = ", "))

if (databases %in% c("stringdb", "ppi", "reactome") && genes_or_proteins == "genes") {
	stop("Can't use a protein based database and map to ensembl gene IDs")
}

# gene names + positions
# ensembl_gene_id hgnc_symbol chromosome_name start_position end_position
all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)
protein_mapping <- read_tsv("data/derived/epigraphdb_uniprot_ids.tsv")

# -------------------------------------------------------
# Setup datasets
# -------------------------------------------------------

all_genes <- all_genes %>%
	mutate(size = end_position - start_position) %>%
	dplyr::filter(!duplicated(ensembl_gene_id))

summary(all_genes$size)

if (genes_or_proteins == "proteins") {
	message("Using protein coding genes mapping")
	all_genes <- all_genes %>%
		dplyr::filter(ensembl_gene_id %in% protein_mapping$ensembl_gene_id) %>%
		dplyr::filter(!duplicated(ensembl_gene_id))
} 


for (db in databases) {
	filenam <- paste0("data/derived/", db, "_terms.txt")
	# read in file
	assign(paste0(db, "_terms"), read_tsv(filenam))
}

# doing this makes performing enrichment analyses much quicker
# databases <- c("kegg", "go", "stringdb")
path <- "data/derived"
db_genes <- lapply(databases, function(db) {
	db_gene_file <- paste0(db, "_genes.rdata")
	full_path <- file.path(path, db_gene_file)
	if (file.exists(full_path)) {
		load(full_path)
		# pathway_genes <- get(paste0(db, "_genes"))
	} else {
		pathway_df <- get(paste0(db, "_terms"))
		pathway_genes <- extract_genes_in_pathways(pathway_df)
		save(pathway_genes, file = full_path)
	}
	# just keep genes in all genes 
	# (important for when using proteins only for pathways that can map to both!) 
	pathway_genes <- map(pathway_genes, function(x) x[x %in% all_genes$ensembl_gene_id])
	return(pathway_genes)
})
names(db_genes) <- databases

for (db in databases) {
	assign(paste0(db, "_genes"), db_genes[[db]])
}

all_unique_genes <- unique(all_genes$ensembl_gene_id)

message("length of all unique genes is ", length(all_unique_genes))

# -------------------------------------------------------
# Setup parameters
# -------------------------------------------------------

# pathway down approach
params1 <- expand.grid(
	prop_causal_pathways = c(0.01),
	prop_consequent_pathways = c(0.005, 0.01, 0.05, 0.1),
	gwas_power = c(0.1),
	ewas_power = c(0.05, 0.1, 0.2, 0.5, 1),
	percent_ewas_causal = c(0, 0.05, 0.1, 0.2, 0.5, 1),
	pathway_type = databases,
	sim = c(1:1000),
	n_gwas_genes = NA,
	n_ewas_genes = NA, 
	gene_overlap = NA, 
	pathway_overlap = NA, 
	or_g = NA,
	or_p = NA, 
	rho_p= NA
)

# gene up approach
params2 <- expand.grid(
	n_genes = c(500, 1000, 2500, 5000, 10000),
	gwas_power = c(0.1),
	ewas_power = c(0.05, 0.1, 0.2, 0.5, 1),
	percent_ewas_causal = c(0, 0.05, 0.1, 0.2, 0.5, 1),
	pathway_type = databases,
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

# params <- params %>%
# 	dplyr::filter(sim == split)

test_params <- params %>%
	dplyr::filter(sim == 1)

# params <- test_params[150:151,]

# have messages set up!
message_at <- calc_increments(nrow(params))

print("starting simulations")

function_name <- paste0("apply_sim", sim_n)
sim_func <- match.fun(function_name)
start_time <- proc.time()
out <- sim_func(params)
out <- bind_rows(out)
time_taken <- proc.time() - start_time
time_taken

path <- "data/derived/sim_temp/"
nam <- paste0("sim", sim_n, "_", paste(databases, collapse = "_"), "_", genes_or_proteins, "_", split, ".txt")
write.table(out, file = paste0(path, nam), 
			row.names = F, col.names = T, quote = F, sep = "\t")


