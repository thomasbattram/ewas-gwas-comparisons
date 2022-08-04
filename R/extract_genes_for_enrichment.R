# ---------------------------------------
# Permuting gene extraction for enrichment tests
# ---------------------------------------

## To assess whether there is more overlap of genes/genesets than 
## expected by chance between GWAS and EWAS of corresponding traits, 
## a null distribution needs to be established using permutations.
## This script extracts genes not related to the trait of interest
## for use in empirical analyses and saves them to save time.

args <- commandArgs(trailingOnly = TRUE)
params_row <- as.numeric(args[1])
traits_to_redo_file <- args[2] ## Only needed if updated EWAS! Doesn't matter about updated GWAS!
home_dir <- args[3] 
# home_dir <- "~/projects/epi_gen_comp"

params <- expand.grid(sampling_method = c("random", "non_random"), 
					  genes_or_proteins = c("genes", "proteins"))

sampling_method <- as.character(params[params_row, "sampling_method"])
genes_or_proteins <- as.character(params[params_row, "genes_or_proteins"])

stopifnot(sampling_method %in% c("random", "non_random"))
stopifnot(genes_or_proteins %in% c("genes", "proteins"))

pkgs <- c("tidyverse", "FDb.InfiniumMethylation.hg19", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "usefunc")
lapply(pkgs, require, character.only = T)

source(file.path(home_dir, "R/mapping_functions.R"))

all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)
cpg_genes <- read_tsv(paste0("data/derived/epigenetic/cpgs_mapped_to_", genes_or_proteins, ".txt"))
snp_genes <- read_tsv(paste0("data/derived/genetic/snps_mapped_to_", genes_or_proteins, ".txt"))
protein_mapping <- read_tsv("data/derived/epigraphdb_uniprot_ids.tsv")

traits <- readLines("data/traits.txt")
traits_to_redo <- readLines(traits_to_redo_file)
# load in data if done previous
out_nam <- paste0("data/derived/", sampling_method, "_permuted_", genes_or_proteins, ".RData")
if (file.exists(out_nam)) {
	load(out_nam)
	old_res <- res
	old_traits <- names(old_res)
	old_traits <- old_traits[!old_traits %in% traits_to_redo]
	traits <- traits[!traits %in% old_traits]
	if (length(traits) == 0) stop("ANALYSIS ALREADY DONE BEFORE!")
}

all_dat <- lapply(traits, function(trait) {
	read_derived_data(trait)
})
names(all_dat) <- traits

if (genes_or_proteins == "proteins") {
	all_genes <- all_genes %>%
		dplyr::filter(ensembl_gene_id %in% protein_mapping$ensembl_gene_id)
}

all_genes <- all_genes %>%
	dplyr::filter(!duplicated(ensembl_gene_id))

max_gene_distance <- max(all_genes$end_position - all_genes$start_position)
min_gene_distance <- min(all_genes$end_position - all_genes$start_position)
chr_positions <- map_df(unique(all_genes$chromosome_name), function(chr) {
	temp <- all_genes %>%
		dplyr::filter(chromosome_name == chr) 
	max_pos <- max(temp$end_position) + 1e6
	min_pos <- min(temp$start_position) - 1e6
	if (min_pos < 1) min_pos <- 1
	out_dat <- data.frame(chromosome = chr, min_pos = min_pos, max_pos = max_pos)
	return(out_dat)
})

# add cumulative position to make it easy to switch between chromosomes
chr_positions <- chr_positions %>%
	rbind(chr_positions) %>%
	rbind(chr_positions) %>%
	rbind(chr_positions) %>%
	dplyr::filter(!chromosome %in% c("X", "Y")) %>%
	mutate(cum_chr = 1:nrow(.)) %>%
	mutate(total_bp = 1 + max_pos - min_pos) %>%
	mutate(cum_pos = cumsum(total_bp)) %>%
	dplyr::select(-total_bp) 

# all_genes %>%
#   mutate(distance = end_position - start_position) %>%
#   dplyr::filter(distance == min(distance))
# checked it out and the length of TRDD1 is actually 7 bases according to genecards...

# input variants linked to genes and output changed genes
extract_genes <- function(variant_dat, permutation, distance, move_chr = TRUE) {
	new_dat <- map_df(unique(variant_dat$chromosome), function(chr) {
	    # new positions for every permutation
		dat <- variant_dat %>%
			dplyr::filter(chromosome == chr) %>%
			mutate(new_position = position + (permutation - 1) * distance)

		# some new positions may be outside the current chromosome bounds
		# so need to move the variants across chromosomes if needs be
		chr_pos <- chr_positions %>%
			dplyr::filter(chromosome == chr) %>%
			dplyr::filter(cum_pos == min(cum_pos))

	    pos_dif <- dat$new_position - chr_pos$max_pos
	    cum_pos_dif <- pos_dif + chr_pos$cum_pos
	    if (any(sign(pos_dif) == 1)) {
		    if (move_chr) {        
		        new_chr <- map_dfr(cum_pos_dif, function(x) {
		        	if (!any(x > chr_positions$cum_pos)) return(data.frame(new_chr = 1, cum_chr = 1))
          			res <- max(which(x > chr_positions$cum_pos)) + 1
          			if (res > 22) {
            			actual_chr <- res - 22
          			} else {
            			actual_chr <- res
          			}
			        # if (res == 23) res <- "X"
			        # if (res == 24) res <- "Y"
			        out_res <- data.frame(new_chr = actual_chr, cum_chr = res)
			        return(out_res)
          			### change this so it outputs the new chromosome and the new cumulative chromosome!!
        		})
		        new_pos <- map_dbl(1:nrow(new_chr), function(x) {
					new_c <- new_chr[x,]
					if (new_c$new_chr == chr) return(dat$new_position[x])          
					new_c_cum_pos <- chr_positions[chr_positions$cum_chr == new_c$cum_chr, "cum_pos", drop = T]
					res <- new_c_cum_pos - cum_pos_dif[x]
					return(res)
		        })
				# change names so that it works with the mapping functions
				dat <- dat %>%
					mutate(position = new_pos) %>%
					mutate(chromosome = as.character(new_chr$new_chr)) %>%
					# mutate(chromosome = case_when(chromosome == "23" ~ "X", 
					#                               chromosome == "24" ~ "Y", 
					#                               !(chromosome %in% c("23", "24")) ~ chromosome)) %>%
					dplyr::select(-new_position)
    			}
    	} else {
      		dat <- dat %>%
        		mutate(position = new_position) %>%
        		dplyr::select(-new_position)
    		}
    	return(dat)
  	})
	# these functions map the new positions to genes
	mapped_dat <- lapply(unique(new_dat$chromosome), map_over_chromosomes, new_dat, all_genes, FALSE)
	mapped_dat <- bind_rows(mapped_dat)
	return(mapped_dat)
}

extract_genes_random <- function(n_variants, permutation, chr_position_dat) {
	set.seed(permutation)
	print(permutation)
	new_cum_pos <- sample(min(chr_position_dat$cum_pos):max(chr_position_dat$cum_pos), 
	               		  n_variants)
	new_dat <- map_dfr(seq_along(new_cum_pos), function(x) {
		temp_cum_pos <- new_cum_pos[x]
		rand_chr <- min(which(temp_cum_pos < chr_position_dat$cum_pos))
		y <- chr_position_dat[rand_chr, "cum_pos"] - temp_cum_pos # how far it is below max
		new_pos <- chr_position_dat[rand_chr, "max_pos"] - y  
		out_dat <- data.frame(name = paste0("cg", x), chromosome = rand_chr, 
		                      position = new_pos)
		return(out_dat)
	})
	print(unique(new_dat$chromosome))
	mapped_dat <- lapply(unique(new_dat$chromosome), map_over_chromosomes, new_dat, all_genes, TRUE)
	mapped_dat <- bind_rows(mapped_dat)
	return(mapped_dat)
}

######

# --------------------------------------------
# Run the extract_genes function on all traits
# --------------------------------------------

start_time <- proc.time()
permutations <- 1001
# permutations <- 6
trait <- traits[1]

old_chr_positions <- chr_positions[1:22, ]
res <- lapply(traits, function(trait) {
	sig_cpgs <- all_dat[[trait]][["meth"]] %>%
		dplyr::filter(P < 1e-7) %>%
		pull(CpG) %>%
    	unique
	
	genes_of_interest <- cpg_genes %>% 
		dplyr::filter(name %in% sig_cpgs)


	out <- lapply(1:permutations, function(perm) {
		print(perm)
	    if (perm == 1) return(genes_of_interest)
	    if (sampling_method == "random") {
	    	print("random sampling of genes!")
	    	GoI <- extract_genes_random(n_variants = length(sig_cpgs), permutation = perm, 
	                                	chr_position_dat = old_chr_positions)
	    } else if (sampling_method == "non_random") {
	    	# extract genes in a way that keeps the correlation structure
	    	GoI <- extract_genes(variant_dat = genes_of_interest, permutation = perm, 
	                             distance = max_gene_distance)
	    }	  	
	    return(GoI)
	})
	return(out)
})
names(res) <- traits

if (exists("old_res")) {
	# out_res <- c(old_res, res)
	res <- c(old_res, res)
} 

fin_time <- proc.time() - start_time
# 71

save(res, file = out_nam)





