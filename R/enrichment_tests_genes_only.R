#
#
# 

all_traits <- c("alcohol_consumption_per_day", 
                "former_versus_never_smoking", 
                "current_versus_never_smoking", 
                "body_mass_index", 
                "glucose", 
                "insulin", 
                "educational_attainment")

databases <- c("kegg", "go", "reactome", "ppi", "stringdb")

sampling_methods <- c("random", "non_random")

params <- expand.grid(trait = all_traits, 
                      database = databases, 
                      sampling_method = sampling_methods)

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = T)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

devtools::load_all("~/repos/usefunc/")

# source in functions for:
# - reading in the data
# - mapping variant positions to genes
source("R/mapping_functions.R")
source("R/simulation_functions.R")

gd <- c("kegg", "go")
# prod <- c("reactome", "ppi", "stringdb")

params$genes_or_proteins <- ifelse(params$database %in% gd, "genes", "proteins")

permuted_proteins_random <- new_load("data/derived/random_permuted_proteins.RData")
permuted_genes_random <- new_load("data/derived/random_permuted_genes.RData")
permuted_proteins_non_random <- new_load("data/derived/non_random_permuted_proteins.RData")
permuted_genes_non_random <- new_load("data/derived/non_random_permuted_genes.RData")

all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)
protein_mapping <- read_tsv("data/derived/epigraphdb_uniprot_ids.tsv")
x=1
start_time <- proc.time()
lapply(1:nrow(params), function(x) {
    message(x)
    df <- params[x, ]
    genes_or_proteins <- df[, "genes_or_proteins", drop=T]
    background_genes <- all_genes
    if (genes_or_proteins == "proteins") {
        background_genes <- background_genes %>%
            dplyr::filter(ensembl_gene_id %in% protein_mapping$ensembl_gene_id)
        snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_proteins.txt")
    } else {
        snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_genes.txt")
    }
    dat <- read_derived_data(df$trait)
    sig_snps <- dat$gen[dat$gen$p < 5e-8, "SNP", drop = T]
    sig_snp_genes <- unique(snp_genes[snp_genes$name %in% sig_snps, "ensembl_gene_id", drop = T])
    rm(snp_genes)
    background_genes <- background_genes %>%
        dplyr::filter(!duplicated(ensembl_gene_id)) %>%
        pull(ensembl_gene_id)
    permuted_genes <- get(paste0("permuted_", df$genes_or_proteins, "_", df$sampling_method))
    permuted_genes <- permuted_genes[[df$trait]]
    # perm <- 1
    message("Starting permutations")
    res <- map_dfr(1:1001, function(perm) {
        GoI <- permuted_genes[[perm]]
        ewas_genes <- unique(GoI$ensembl_gene_id)
        gene_res <- overlap_test(group1 = ewas_genes, 
                                 group2 = sig_snp_genes, 
                                 all_variables = background_genes)
        gene_out_res <- tibble(perm = perm, OR = gene_res$estimate, p = gene_res$p)
        return(gene_out_res)
    })
    pre_nam <- paste0("results/", df$trait, "/", df$trait, "_", df$sampling_method, "_sampling_", df$database, "_")
    nam <- paste0(pre_nam, "gene_enrichment_tests.txt")
    message("writing out table")
    write.table(res, file = nam, 
                row.names = F, col.names = T, quote = F, sep = "\t")
    message("FIN!")
})
time_taken <- proc.time() - start_time
time_taken
