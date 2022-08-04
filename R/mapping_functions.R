# ------------------------------------------------------------
# file with functions in used in multiple scripts for this project
# ------------------------------------------------------------

#' Map cpg positions or snp positions to genes
#' 
#' @param dat data.frame with the variants of interest and their position in it
#' @param gene_dat data.frame with all the gene start and end positions in
#' @return a data.frame with the variants, their positions, and the genes they map to (ensembl id and hgnc)
map_positions_to_genes <- function(dat, gene_dat) 
{
	# sort the CpGs within genes first
	new_gene_dat <- gene_dat %>%
		dplyr::filter(!end_position < min(dat$position)) %>%
		dplyr::filter(!start_position > max(dat$position))

	if (nrow(new_gene_dat) == 0) return(NULL)

	mapped_variants <- lapply(1:nrow(new_gene_dat), function(x) {
		gene <- new_gene_dat[x, ]
		variants_within_gene <- between(dat$position, gene$start_position, gene$end_position)
		if (sum(variants_within_gene) == 0) return(NULL)		
		variant_dat <- dat %>%
			dplyr::filter(variants_within_gene) %>%
			mutate(ensembl_gene_id = gene$ensembl_gene_id, hgnc_symbol = gene$hgnc_symbol,
				   distance_to_gene = 0) %>%
			distinct()
		return(variant_dat)
	})
	within_gene_dat <- bind_rows(mapped_variants)
	
	# then map the rest to the nearest gene! 
	unmapped_variants <- dat %>%
		dplyr::filter(!name %in% within_gene_dat$name)
	if (nrow(unmapped_variants) == 0) return (within_gene_dat)
	
	new_gene_dat <- gene_dat %>%
		gather(key = where, value = position, -one_of(c("ensembl_gene_id", "hgnc_symbol", "uniprot_gn_id", "chromosome_name")))

	outside_gene_dat <- lapply(1:nrow(unmapped_variants), function(x) {
		temp <- unmapped_variants[x, ]
		gene_distances <- abs(new_gene_dat$position - abs(temp$position))
		closest_gene_distance <- min(gene_distances)
		closest_gene <- which.min(gene_distances)
		temp_gen <- new_gene_dat[closest_gene,]
		out_dat <- temp %>%
			mutate(ensembl_gene_id = temp_gen$ensembl_gene_id, hgnc_symbol = temp_gen$hgnc_symbol, 
				   distance_to_gene = closest_gene_distance)
		return(out_dat)
	})
	outside_gene_dat <- bind_rows(outside_gene_dat)

	# bind the two together and output it! 
	out_dat <- rbind(within_gene_dat, outside_gene_dat)
	return(out_dat)
}

#' Map positions to genes over chromosomes
#' 
#' @param chr chromosome number
#' @param variant_dat data.frame with variant names and positions in 
#' @param gene_dat data.frame with all genes with their start and end positions
#' @return list with results from map_positions_to_genes function across all chromosomes
map_over_chromosomes <- function(chr, variant_dat, gene_dat, verbose = TRUE) 
{
	if (verbose) print(chr)
	if (!"chromosome" %in% colnames(variant_dat)) {
		chr_col_name <- grep("^chr", colnames(variant_dat), ignore.case = TRUE, value = TRUE)
		stopifnot(length(chr_col_name) == 1)
		colnames(variant_dat)[colnames(variant_dat) == chr_col_name] <- "chromosome"
	}
	dat <- dplyr::filter(variant_dat, chromosome == chr)
	gene_dat <- dplyr::filter(gene_dat, chromosome_name == chr)
	fin <- map_positions_to_genes(dat, gene_dat)
	return(fin)
}

#' Read in sorted genetic and DNA methylation data from file for a given trait
#' 
#' @param trait the trait name which should also be in the name of the file
#' @return a list containing the sorted genetic and DNA methylation data
read_derived_data <- function(trait) 
{
	file_path <- "data/derived/"
	gen_dat <- read_tsv(paste0(file_path, "genetic/sorted_", trait, ".txt"), guess_max = 1e6)
	meth_dat <- read_tsv(paste0(file_path, "epigenetic/sorted_", trait, ".txt"), guess_max = 1e6)
	return(list(gen = gen_dat, meth = meth_dat))
}

#' Check the results directories are there
#' 
#' @return TRUE or FALSE depending on directory is there
check_dir <- function(trait) 
{
	file_path <- "results/"
	files <- list.files(file_path)
	out <- TRUE
	if (!trait %in% files) {
		message(paste0("There is no ", trait, " results directory"))
		out <- FALSE
	}
	return(out)
}

#' Checks all the results directories are there
#' 
#' @return TRUE or FALSE depending on whether the directories are there
check_all_dirs <- function(traits) 
{
	out <- purrr::map_lgl(traits, function(trait) {
		check_dir(trait)
	})
	stopifnot(all(out))
}

