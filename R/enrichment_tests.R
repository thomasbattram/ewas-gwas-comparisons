# ------------------------------------------------------------------------
# Testing overlap of EWAS/GWAS genes + pathways
# ------------------------------------------------------------------------

### Script aims to test the overlap between genes and pathways identified
### by corresponding GWAS and EWAS using Fisher's exact test.

args <- commandArgs(trailingOnly = TRUE)
params_row <- as.numeric(args[1])
traits_to_redo_file <- args[2]

home_dir <- args[3]
# home_dir <- "~/projects/epi_gen_comp" # CHANGE ME WHEN NEEDED

all_traits <- readLines("data/traits.txt")
traits_to_redo <- readLines(traits_to_redo_file)
databases <- c("go") ## can add more in if want to, but only go used for paper!

sampling_methods <- c("random", "non_random")

params <- expand.grid(trait = all_traits, 
                      database = databases, 
                      sampling_method = sampling_methods)

# params$gene_enrich <- ifelse(params$database == "kegg" | params$database == "reactome", TRUE, FALSE)

trait <- as.character(params[params_row, "trait"])
path_type <- as.character(params[params_row, "database"])
sampling_method <- as.character(params[params_row, "sampling_method"])
# gene_enrich <- params[params_row, "gene_enrich"]

message("The trait is ", trait)
message("The pathway is ", path_type)
message("The sampling method is ", sampling_method)

stopifnot(path_type %in% c("kegg", "go", "reactome", "ppi", "stringdb"))
stopifnot(sampling_method %in% c("random", "non_random"))

library(tidyverse) # tidy code and data
library(usefunc) # personal package of useful functions

# source in functions for:
# - reading in the data
# - mapping variant positions to genes
source(file.path(home_dir, "R/mapping_functions.R"))
source(file.path(home_dir, "R/simulation_functions.R"))

check_all_dirs(trait)

#
# test with just BMI data
#

# ------------------------------------------------------------------------
# If analysis not done before then read in data
# ------------------------------------------------------------------------

# check if all output files expected are present
pre_nam <- paste0("results/", trait, "/", trait, "_", sampling_method, "_sampling_", path_type, "_")

nams <- c("pathway_enrichment_tests.txt", "pathway_enrichment_score_correlations.txt", 
          "pathway_and_gene_overlap.txt", "genes_from_pathway_overlap_tests.RData", 
          "pathways_from_pathway_overlap_tests.RData")

if (!trait %in% traits_to_redo) {
    lapply(nams, function(x) {
        file_path <- paste0(pre_nam, x)
        if (file.exists(file_path)) stop("Quiting because files already exist!")
    })
}
# read in data
dat <- read_derived_data(trait)
pathway_dat_nam <- paste0("data/derived/", path_type, "_terms.txt")
pathway_dat <- read_tsv(pathway_dat_nam)
all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)
#
#
#

if (path_type %in% c("kegg", "go")) {
    cpg_genes <- read_tsv("data/derived/epigenetic/cpgs_mapped_to_genes.txt")
    snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_genes.txt")
    all_genes <- all_genes %>%
        dplyr::filter(!duplicated(ensembl_gene_id))
    permuted_genes_nam <- paste0("data/derived/", sampling_method, "_permuted_genes.RData")
    load(permuted_genes_nam)
} else {
    cpg_genes <- read_tsv("data/derived/epigenetic/cpgs_mapped_to_proteins.txt")
    snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_proteins.txt")
    protein_mapping <- read_tsv("data/derived/epigraphdb_uniprot_ids.tsv")
    all_genes <- all_genes %>%
        dplyr::filter(ensembl_gene_id %in% protein_mapping$ensembl_gene_id) %>%
        dplyr::filter(!duplicated(ensembl_gene_id))
    permuted_genes_nam <- paste0("data/derived/", sampling_method, "_permuted_proteins.RData")
    load(permuted_genes_nam)
}
permuted_genes <- res[[trait]]
rm(res)

# ---------------------------------------------
# Setup analysis
# ---------------------------------------------

sig_cpgs <- dat$meth[dat$meth$P < 1e-7, "CpG", drop = T]
genes_of_interest <- cpg_genes[cpg_genes$name %in% sig_cpgs, ]
# genes_of_interest <- cpg_genes[cpg_genes$name %in% sig_cpgs, "ensembl_gene_id", drop = T]
sig_snps <- dat$gen[dat$gen$p < 5e-8, "SNP", drop = T]
sig_snp_genes <- unique(snp_genes[snp_genes$name %in% sig_snps, "ensembl_gene_id", drop = T])

# extract only the go terms needed --> for speed! 
new_pathway_dat <- pathway_dat %>%
    dplyr::filter(!is.na(pathway_id))

# remove unnecessary large clutter
rm(snp_genes)

# this makes enrichment much quicker!
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


pathways <- unique(new_pathway_dat$pathway_id)

gwas_pathways <- new_pathway_dat %>%
    dplyr::filter(ensembl_gene_id %in% sig_snp_genes) %>%
    pull(pathway_id) %>%
    unique

# ------------------------------------------------------------------------
# Run analysis
# ------------------------------------------------------------------------

correlate_enrichment_scores <- function(gwas_scores, ewasg, cor_test = FALSE)
{
    ### function to correlate the enrichment scores
    ewas_enrich <- perform_enrichment(identified_genes = ewasg, 
                                      database = path_type, 
                                      background_genes = all_genes$ensembl_gene_id)
    
    gwas_s <- log(gwas_scores$enrich_or)
    ewas_s <- log(ewas_enrich$enrich_or)
    if (cor_test) {
        enrich_cor_test <- cor.test(gwas_s, ewas_s, use = "pairwise.complete.obs", 
                               method = "spearman")
        enrich_cor <- enrich_cor_test$estimate
        if (is.na(enrich_cor)) {
            enrich_cor <- 0
            names(enrich_cor) <- 1
        } else {
            names(enrich_cor) <- comma(enrich_cor_test$p.value)
        }
    } else {
        enrich_cor <- cor(gwas_s, ewas_s,
                          use = "pairwise.complete.obs",
                          method = "spearman")        
    }
    return(enrich_cor)
}

gwas_enrich <- perform_enrichment(sig_snp_genes, 
                                  database = path_type, 
                                  background_genes = all_genes$ensembl_gene_id)

ewas_genes <- unique(permuted_genes[[1]]$ensembl_gene_id)
enrich_cor <- correlate_enrichment_scores(gwas_scores=gwas_enrich,
                                          ewasg = ewas_genes, 
                                          cor_test=TRUE)
# Cannot compute exact p-value with ties error will come up because
# there are plenty of times where the enrichment score will be zero
# for both gwas and ewas.

# ewasee <- arrange(ewas_enrich, desc(enrich_or)) %>%
#     dplyr::filter(enrich_or != Inf & enrich_or != 0) %>%
#     as_tibble
# gwasee <- arrange(gwas_enrich, desc(enrich_or)) %>%
#     dplyr::filter(enrich_or != Inf & enrich_or != 0) %>%
#     dplyr::filter(pathway_id %in% ewasee$pathway_id) %>%
#     as_tibble 
# ewasee <- dplyr::filter(ewasee, pathway_id %in% gwasee$pathway_id)


# cor(log(gwasee$enrich_or), log(ewasee$enrich_or),
#                           use = "pairwise.complete.obs",
#                           method = "spearman") 



cur_time <- proc.time()
# perm=1
res <- lapply(1:1001, function(perm) {
    print(perm)
    GoI <- permuted_genes[[perm]]
    ewas_genes <- unique(GoI$ensembl_gene_id)
    # map GoIs to pathways 
    PoI <- new_pathway_dat %>%
        dplyr::filter(ensembl_gene_id %in% ewas_genes)

    ewas_pathways <- unique(PoI[["pathway_id"]])
    ewas_genes_in_pathways <- unique(PoI[["ensembl_gene_id"]])
    n_pathway_overlap <- sum(ewas_pathways %in% gwas_pathways)

    # perform the overlap
    path_res <- overlap_test(group1 = ewas_pathways, 
                             group2 = gwas_pathways, 
                             all_variables = pathways)

    # perform enrichment score tests
    if (as.numeric(names(enrich_cor)) < 0.05) {
        cor_res <- correlate_enrichment_scores(gwas_enrich, ewas_genes)
        enrich_out_res <- data.frame(perm = perm, cor = cor_res)
    } else {
        enrich_out_res <- NULL
    }

    path_out_res <- data.frame(perm = perm, OR = path_res$estimate, p = path_res$p)
    # now for genes
    gene_res <- overlap_test(group1 = ewas_genes, 
                             group2 = sig_snp_genes, 
                             all_variables = all_genes$ensembl_gene_id)

    gene_out_res <- data.frame(perm = perm, OR = gene_res$estimate, p = gene_res$p)

    n_gene_overlap <- sum(ewas_genes %in% sig_snp_genes)

    # see how many pathways overlap with the genetic pathways identified
    overlap_res <- data.frame(perm = perm, pathway_overlap = n_pathway_overlap, n_ewas_pathways = length(ewas_pathways), 
                            n_ewas_genes = length(unique(ewas_genes)), gene_overlap = n_gene_overlap) %>%
                   mutate(prop_pathway_overlap = pathway_overlap / n_ewas_pathways, genes_in_pathways = length(ewas_genes_in_pathways))

    out_res <- list(pathway_res = path_out_res,
                    gene_res = gene_out_res, 
                    enrich_res = enrich_out_res,
                    overlap = overlap_res, 
                    genes = GoI, 
                    pathways = PoI)
    return(out_res)
})

pathway_res <- map_dfr(res, "pathway_res")
overlap_res <- map_dfr(res, "overlap")
enrich_res <- map_dfr(res, "enrich_res")
gene_dat <- map(res, "genes")
pathway_dat <- map(res, "pathways")
gene_res <- map_dfr(res, "gene_res")

fin_time <- proc.time() - cur_time
fin_time 
# KEGG: 3.9s, thus ~65mins for all 
# GO: 93s, thus 26h for all
# reactome: 5.6s, thus ~93mins for all
# ppi: 78s, thus 22h for all
pre_nam <- paste0("results/", trait, "/", trait, "_", sampling_method, "_sampling_", path_type, "_")

nam <- paste0(pre_nam, "pathway_enrichment_tests.txt") 
write.table(pathway_res, file = nam, 
            row.names = F, col.names = T, quote = F, sep = "\t")

nam <- paste0(pre_nam, "gene_enrichment_tests.txt")
write.table(gene_res, file = nam, 
            row.names = F, col.names = T, quote = F, sep = "\t")

if (nrow(enrich_res) == 0) {
    enrich_res <- data.frame(cor_est = enrich_cor, cor_p = names(enrich_cor))
}

nam <- paste0(pre_nam, "pathway_enrichment_score_correlations.txt")
write.table(enrich_res, file = nam, 
            row.names = F, col.names = T, quote = F, sep = "\t")

nam <- paste0(pre_nam, "pathway_and_gene_overlap.txt")
write.table(overlap_res, file = nam, 
            row.names = F, col.names = T, quote = F, sep = "\t")

nam <- paste0(pre_nam, "genes_from_pathway_overlap_tests.RData")
save(gene_dat, file = nam)

nam <- paste0(pre_nam, "pathways_from_pathway_overlap_tests.RData")
save(pathway_dat, file = nam)

# renaming files
# traits <- c("alcohol_consumption_per_day", "former_versus_never_smoking", "current_versus_never_smoking",
#      "body_mass_index", "glucose", "insulin", "educational_attainment")
# pathway_db <- c("gene_ontology", "kegg")
# sampling_method <- "non_random"
# tra=traits[1]
# path_type=pathway_db[1]
# lapply(traits, function(tra) {
#   print(tra)
#   path <- paste0("results/", tra, "/")
#   file_nam <- list.files(path)
#   lapply(pathway_db, function(path_type) {
#     old_pre_nam <- paste0(tra, "_" , path_type, "_")
#     pre_nam <- paste0(tra, "_", sampling_method, "_sampling_", path_type, "_")
#     old_nam <- grep(old_pre_nam, file_nam, value = TRUE)
#     new_nam <- gsub(old_pre_nam, pre_nam, old_nam)
    
#     lapply(seq_along(old_nam), function(x) {
#       cmd <- paste0("mv ", path, old_nam[x], " ", path, new_nam[x])
#       system(cmd)
#       return(NULL)
#     })
#     return(NULL)
#   })
#   return(NULL)
# })


