# -------------------------------------------------------
# Simulating EWAS and GWAS pathway overlap
# -------------------------------------------------------

## Simulations that take empirically identified genes/genesets and
## ask if these can be used to assess likely overlap of causal and associated
## genes given varying levels of genes to discover

pkgs <- c("tidyverse", "FDb.InfiniumMethylation.hg19",
		  "IlluminaHumanMethylation450kanno.ilmn12.hg19", "usefunc")

lapply(pkgs, require, character.only = T)

## args
args <- commandArgs(trailingOnly = TRUE)
# args <- c('100', 'egfr')
split <- as.numeric(args[1])
split1 <- (split - 1) * 10 + 1
split2 <- split * 10
message("split = ", split1, " to ", split2)
trait <- args[2]
home_dir <- args[3]
# home_dir <- "~/projects/epi_gen_comp/"

source(file.path(home_dir, "R/mapping_functions.R"))
source(file.path(home_dir, "R/simulation_functions.R"))
# gene ontology terms - table of terms + genes linked to the terms
# ensembl_gene_id pathway_id description
pathway_db <- "go"

go_terms <- read_tsv("data/derived/go_terms.txt") %>% 
	dplyr::filter(!is.na(pathway_id))

# gene names + positions
# ensembl_gene_id hgnc_symbol chromosome_name start_position end_position
all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)
all_genes$size <- all_genes$end_position - all_genes$start_position
summary(all_genes$size)

# unique pathways
u_go <- unique(go_terms$pathway_id)

# -------------------------------------------------------
# Extract extra information for parameters 3 & 4
# -------------------------------------------------------
traits <- readLines("data/traits.txt")
if (!trait %in% traits) stop("Trait entered is faulty")
pathway_databases <- c("go")
sampling_method <- "random"

# extract ewas genes here
ewas_genes_file <- "data/derived/ewas_genes.txt"
if (file.exists(ewas_genes_file)) {
	ewas_g_dat <- read_tsv(ewas_genes_file)
}

if (any(!traits %in% ewas_g_dat$trait) | !file.exists(ewas_genes_file)) {
	cpg_genes <- read_tsv("data/derived/epigenetic/cpgs_mapped_to_genes.txt")
	ewas_g_dat <- map_dfr(traits, function(trait) {
		dat <- read_derived_data(trait)
		sig_cpgs <- dat$meth[dat$meth$P < 5e-8, "CpG", drop = T]
		sig_cpg_genes <- unique(cpg_genes[cpg_genes$name %in% sig_cpgs, "ensembl_gene_id", drop = T])
		out <- data.frame(trait = trait, ewas_genes = unique(sig_cpg_genes))
		return(out)
	})
	write.table(ewas_g_dat, file = ewas_genes_file, 
			   row.names = F, col.names = T, quote = F, sep = "\t")
	rm(cpg_genes)
}

# extract gwas genes here 
gwas_genes_file <- "data/derived/gwas_genes.txt"
if (file.exists(gwas_genes_file)) {
	gwas_g_dat <- read_tsv(gwas_genes_file)
}

if (any(!traits %in% gwas_g_dat$trait) | !file.exists(gwas_genes_file)) {
	snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_genes.txt")
	gwas_g_dat <- map_dfr(traits, function(trait) {
		dat <- read_derived_data(trait)
		sig_snps <- dat$gen[dat$gen$p < 5e-8, "SNP", drop = T]
		sig_snp_genes <- unique(snp_genes[snp_genes$name %in% sig_snps, "ensembl_gene_id", drop = T])
		out <- data.frame(trait = trait, gwas_genes = unique(sig_snp_genes))
		return(out)
	})
	write.table(gwas_g_dat, file = gwas_genes_file, 
			   row.names = F, col.names = T, quote = F, sep = "\t")
	rm(snp_genes)
}

trait_gene_dat <- ewas_g_dat %>%
	left_join(gwas_g_dat) %>%
	dplyr::select(trait, ewas_genes, gwas_genes) %>%
	distinct() 

trait_n_gene_dat <- trait_gene_dat %>%
	group_by(trait) %>%
	summarise(n_ewas_genes = length(unique(ewas_genes)), n_gwas_genes = length(unique(gwas_genes)))

trait_n_gene_dat %>%
	mutate(prop_ewas = n_ewas_genes / nrow(all_genes), prop_gwas = n_gwas_genes / nrow(all_genes))

# minimum prop_causal_genes or prop_consequential_genes = 0.04 
# (any lower would mean the number of genes identified is greater than present in sample)

# get genes to make things quicker for running enrichment tests!
path <- "data/derived"
db_gene_file <- paste0(pathway_databases, "_genes.rdata")
full_path <- file.path(path, db_gene_file)
if (file.exists(full_path)) {
    load(full_path)
    assign(paste0(pathway_databases, "_genes"), pathway_genes)
} else {
    pathway_df <- get(paste0(db, "_terms"))
    pathway_genes <- extract_genes_in_pathways(pathway_df)
    save(pathway_genes, file = full_path)
    assign(paste0(pathway_databases, "_genes"), pathway_genes)
}


# -------------------------------------------------------
# Setup parameters for simulations 3 & 4
# -------------------------------------------------------

params <- expand.grid(
	trait = trait,
	n_gene_proportions =  c(1, 2, 3, 5, 10, 20), 
	ca_con_overlap = c(0, 0.01, 0.1, 0.5, 1),
	sim = c(1:1000),
	gene_overlap = NA, 
	or_g = NA, 
	or_p = NA, 
	rho_p = NA
	)

# -------------------------------------------------------
# run the simulations
# -------------------------------------------------------

params <- params %>%
	dplyr::filter(between(sim, split1, split2))

test_params <- params[37:38,]

# params <- test_params

print("starting simulations")

start_time <- proc.time()
out <- apply_sim4(params)
out <- bind_rows(out)
time_taken <- proc.time() - start_time
time_taken

path <- "data/derived/sim_temp/architecture_sims/"
nam <- paste0("architecture_sims_", pathway_db, "_", trait, "_", split, ".txt")
write.table(out, file = paste0(path, nam), 
			row.names = F, col.names = T, quote = F, sep = "\t")


