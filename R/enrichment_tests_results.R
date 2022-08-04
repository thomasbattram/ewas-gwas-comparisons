# ------------------------------------------------------------------
# Taking the results from the pathway analysis and visualising them
# ------------------------------------------------------------------

## Aims of script:
## 1. Explore enrichment results
## 2. Synthesise results in a way that can be presented in a paper

## pkgs
library(tidyverse) # tidy code + data
library(gridExtra) # plots on multiple pages -- marrangeGrob()
library(usefunc) # personal package of useful functions

# setwd("SCRATCH_SPACE")
home_dir <- "~/projects/epi_gen_comp" # CHANGE ME WHEN NEEDED

source(file.path(home_dir, "R/mapping_functions.R"))

## data
traits <- readLines("data/traits.txt")
pathway_databases <- c("go")
sampling_methods <- c("random", "non_random")

## bit messy script, so here is structure:
## 1. Comparison of random and non-random sampling of positions in genome
##	  -- this section essentially compares two methods of generating a null
##	     distribution. One using random genes and another aiming to keep the
##		 "structure" in DNA methylation data
## 2. Is there more pathway correlation when using the original EWAS results than by chance
##	  -- this section checks the distribution of permuted results and compares the permuted
##		 results (null distribution) to the empirical results. It also looks at whether the
##	  	 number of pathways may contribute to the results
## 3. Is there more gene overlap when using the original EWAS results than by chance
##	  -- same as 2, but with gene overlap results
## 4. What pathways tend to come up all the time?
## 	  -- Checks whether there are certain pathways that are identified regularly across
##		 permutations and traits
## 5. Extraction of data for paper
##	  -- generation of tables for the paper (bit messy)


# ----------------------------------------------------
# comparison of random and non-random sampling of positions in genome
# ----------------------------------------------------

temp <- expand.grid(trait = traits,
 					pd = pathway_databases,
 					sampling_method = sampling_methods)

trait_res <- map_dfr(1:nrow(temp), function(row) {
	df <- temp[row, ]
	print(row)
	pre_nam <- paste0("results/", df$trait, "/",
					  df$trait, "_", df$sampling_method, "_sampling_",  df$pd, "_")
	pathway_nam <- paste0(pre_nam, "pathway_enrichment_tests.txt")
	pathway_out <- suppressMessages(read_tsv(pathway_nam)) %>%
		mutate(trait = df$trait, pd = df$pd, sm = df$sampling_method) %>%
		dplyr::select(trait, pd, sm, everything())
	enrich_nam <- paste0(pre_nam, "pathway_enrichment_score_correlations.txt")
	enrich_out <- suppressMessages(read_tsv(enrich_nam)) %>%
		mutate(trait = df$trait, pd = df$pd, sm = df$sampling_method) %>%
		dplyr::select(trait, pd, sm, everything())	
	overlap_nam <- paste0(pre_nam, "pathway_and_gene_overlap.txt")
	overlap_out <- suppressMessages(read_tsv(overlap_nam)) %>%
		mutate(trait = df$trait, pd = df$pd, sm = df$sampling_method) %>%
		dplyr::select(trait, pd, sm, everything())
	comb_res <- pathway_out %>%
		left_join(overlap_out) %>%
		left_join(enrich_out)
	return(comb_res)
})

perm_only_res <- trait_res %>%
	dplyr::filter(perm != 1)

# x = traits
# y = value
# fill = sampling_method
# just use gene ontology

vals <- c("OR", "pathway_overlap", "gene_overlap")

p_list <- lapply(seq_along(vals), function(x) {
	val <- vals[x]
	p_dat <- perm_only_res %>%
		dplyr::filter(pd == "gene_ontology") %>%
		dplyr::select(trait, sm, one_of(val))
	p_out <- ggplot(p_dat, aes_string(x = "trait", y = val, fill = "sm")) +
		geom_boxplot()
	return(p_out)
})

glist <- lapply(p_list, ggplotGrob)
ggsave("results/plots/comparison_of_position_sampling_methods.pdf", 
		plot = marrangeGrob(glist, nrow=1, ncol=1, top = NULL))


# ------------------------------------------------------------------
# is there more pathway correlation when using the original EWAS results than by chance
# ------------------------------------------------------------------
temp <- expand.grid(trait = traits,
 					pd = pathway_databases,
 					sampling_method = sampling_methods)

res <- map_dfr(1:nrow(temp), function(row) {
	df <- temp[row, ]
	nam <- paste0("results/", df$trait, "/", df$trait, "_", df$sampling_method, "_sampling_", df$pd, "_", "pathway_enrichment_score_correlations.txt")
	out <- read_tsv(nam) %>%
		mutate(trait = df$trait, pd = df$pd, sm = df$sampling_method) %>%
		dplyr::select(trait, pd, sm, everything())
	return(out)
})

perm1_res <- res %>%
	dplyr::filter(perm == 1)

# looking at the data
res %>%
	group_by(trait, pd, sm) %>%
	summarise(med_cor = median(cor), max_cor = max(cor)) %>%
	left_join(perm1_res) %>%
	dplyr::select(-perm)
### so essentially the real results are no different to the null...

# what does the correlation distribution look like? -- checking null dist!
p <- lapply(1:nrow(temp), function(row) {
	temp_trait <- temp[row, "trait"]
	temp_pd <- temp[row, "pd"] 
	df <- res %>%
		dplyr::filter(trait == temp_trait & pd == temp_pd)
	out <- ggplot(df, aes(x = cor)) +
		geom_histogram() +
		theme_bw() +
		labs(title = paste0(temp_trait, "_", temp_pd))
	return(out)
})

glist <- lapply(p, ggplotGrob)
ggsave("results/plots/pathway_overlap_OR_distribution.pdf", 
		plot = marrangeGrob(glist, nrow=1, ncol=1, top = NULL))


temp <- expand.grid(trait = traits,
 					pd = pathway_databases,
 					sampling_method = sampling_methods)
overlap_res <- map_dfr(1:nrow(temp), function(row) {
	df <- temp[row, ]
	nam <- paste0("results/", df$trait, "/", df$trait, "_", df$sampling_method, "_sampling_", df$pd, "_", "pathway_and_gene_overlap.txt")
	out <- read_tsv(nam) %>%
		mutate(trait = df$trait, pd = df$pd) %>%
		dplyr::select(trait, pd, everything())
	return(out)
})


# Does cor also associate with number of EWAS pathways?
comb_res <- res %>%	
	left_join(overlap_res)

pathway_p_list <- lapply(1:nrow(temp), function(row) {
	df <- temp[row, ]
	p_res <- comb_res %>%
		dplyr::filter(trait == df$trait & pd == df$pd)
	out1 <- ggplot(p_res, aes(x = perm, y = n_ewas_pathways)) +
		geom_point() +
		theme_bw() +
		labs(title = paste0(df$trait, "_", df$pd))

	out2 <- ggplot(p_res, aes(x = n_ewas_pathways, y = cor)) +
		geom_point() +
		theme_bw()
	return(list(out1, out2))
})

p_list2 <- flatten(pathway_p_list)

glist <- lapply(p_list2, ggplotGrob)
ggsave("results/plots/pathway_enrichment_correlation_and_n_ewas.pdf",
	   plot = marrangeGrob(glist, nrow = 1, ncol = 2, top = NULL))

# ------------------------------------------------------------------
# is there more gene overlap when using the original EWAS results than by chance
# ------------------------------------------------------------------

temp <- expand.grid(trait = traits,
 					pd = pathway_databases,
 					sampling_method = sampling_methods) 
res <- map_dfr(1:nrow(temp), function(row) {
	df <- temp[row, ]
	nam <- paste0("results/", df$trait, "/", df$trait, "_", df$sampling_method, "_sampling_", df$pd, "_gene_enrichment_tests.txt")
	out <- read_tsv(nam) %>%
		mutate(trait = df$trait) %>%
		dplyr::select(trait, everything())
	return(out)
})

perm1_res <- res %>%
	dplyr::filter(perm == 1)

# looking at the data
res %>%
	group_by(trait) %>%
	summarise(med_p = median(p), med_or = median(OR), max_or = max(OR)) %>%
	left_join(perm1_res) %>%
	dplyr::select(-perm)
### so essentially the real results are no different to the null...

# what does the OR value distribution look like?
p <- lapply(1:nrow(temp), function(row) {
	temp_trait <- temp[row, "trait"]
	df <- res %>%
		dplyr::filter(trait == temp_trait)
	out <- ggplot(df, aes(x = OR)) +
		geom_histogram() +
		theme_bw() +
		labs(title = temp_trait)
	return(out)
})

glist <- lapply(p, ggplotGrob)
ggsave("results/plots/gene_overlap_OR_distribution.pdf", 
		plot = marrangeGrob(glist, nrow=1, ncol=1, top = NULL))

temp <- expand.grid(trait = traits,
 					pd = pathway_databases,
 					sampling_method = sampling_methods)
overlap_res <- map_dfr(1:nrow(temp), function(row) {
	df <- temp[row, ]
	nam <- paste0("results/", df$trait, "/", df$trait, "_", df$sampling_method, "_sampling_", df$pd, "_", "pathway_and_gene_overlap.txt")
	out <- read_tsv(nam) %>%
		mutate(trait = df$trait, pd = df$pd) %>%
		dplyr::select(trait, pd, everything())
	return(out)
})

# Does OR also associate with number of EWAS pathways?
comb_res <- res %>%	
	left_join(overlap_res)

p_list <- lapply(1:nrow(temp), function(row) {
	df <- temp[row, ]
	p_res <- comb_res %>%
		dplyr::filter(trait == df$trait & pd == df$pd)
	out1 <- ggplot(p_res, aes(x = perm, y = n_ewas_pathways)) +
		geom_point() +
		theme_bw() +
		labs(title = paste0(df$trait, "_", df$pd))

	out2 <- ggplot(p_res, aes(x = n_ewas_pathways, y = OR)) +
		geom_point() +
		theme_bw()
	return(list(out1, out2))
})

p_list2 <- flatten(p_list)

glist <- lapply(p_list2, ggplotGrob)
ggsave("results/plots/overlap_and_n_ewas.pdf",
	   plot = marrangeGrob(glist, nrow = 1, ncol = 2, top = NULL))

# -------------------------------------------------------
# what pathways tend to come up all the time?
# -------------------------------------------------------

temp <- expand.grid(trait = traits,
 					pd = pathway_databases,
 					sampling_method = sampling_methods)
x=1
pathway_res <- map_dfr(1:nrow(temp), function(x) {
	print(x)
	df <- temp[x, ]
	nam <- paste0("results/", df$trait, "/", df$trait, "_", df$sampling_method, "_sampling_", df$pd, "_", "pathways_from_pathway_overlap_tests.RData")
	load(nam)
	pathway_res <- map_dfr(seq_along(pathway_dat), function(x) {
		dat <- pathway_dat[[x]] %>%
			dplyr::filter(!duplicated(pathway_id)) %>%
			dplyr::select(-ensembl_gene_id)
		return(dat)
	})	
	
	pathway_res2 <- pathway_res %>%
		group_by(pathway_id) %>%
		summarise(count = n()) %>%
		arrange(desc(count)) %>%
		dplyr::filter(count > 500) %>%
		# head(n = 10) %>%
		mutate(trait = df$trait, pd = df$pd)

	return(pathway_res2)
})

sum_p_res <- pathway_res %>%
	group_by(pathway_id) %>%
	summarise(total_count = sum(count), traits = n()) %>%
	arrange(desc(total_count))


summary(sum_p_res)
# most of the time if a pathway is identified lots for permutations
# for one trait, it doesn't tend to be found for the others...
sum_p_res %>%
	dplyr::filter(traits == 15)
# only 72 out of 3588 pathways are regularly present in over half of 
# permutations for all traits 

# would expect normal distribution, but not case for all of them... 

# ------------------------------------------------------
# Extraction of data for paper
# ------------------------------------------------------
# getting n_gwas_genes
n_gwas_genes_file <- "data/derived/n_gwas_genes.txt"
if (file.exists(n_gwas_genes_file)) {
	gwas_g_dat <- read_tsv(n_gwas_genes_file)
	if (!all(traits %in% gwas_g_dat$trait)) warning("NOT ALL TRAITS PRESENT")
} else {
	snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_genes.txt")
	gwas_g_dat <- map_dfr(traits, function(trait) {
		dat <- read_derived_data(trait)
		sig_snps <- dat$gen[dat$gen$p < 5e-8, "SNP", drop = T]
		sig_snp_genes <- unique(snp_genes[snp_genes$name %in% sig_snps, "ensembl_gene_id", drop = T])
		out <- data.frame(trait = trait, n_gwas_genes = length(unique(sig_snp_genes)))
		return(out)
	})
	write.table(gwas_g_dat, file = n_gwas_genes_file, 
			   row.names = F, col.names = T, quote = F, sep = "\t")
	rm(snp_genes)
}

# getting n_gwas_proteins 
n_gwas_proteins_file <- "data/derived/n_gwas_proteins.txt"
if (file.exists(n_gwas_proteins_file)) {
	gwas_p_dat <- read_tsv(n_gwas_proteins_file)
	if (!all(gwas_g_dat$trait %in% traits)) warning("NOT ALL TRAITS PRESENT")
} else {
	snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_proteins.txt")
	gwas_p_dat <- map_dfr(traits, function(trait) {
		dat <- read_derived_data(trait)
		sig_snps <- dat$gen[dat$gen$p < 5e-8, "SNP", drop = T]
		sig_snp_genes <- unique(snp_genes[snp_genes$name %in% sig_snps, "ensembl_gene_id", drop = T])
		out <- data.frame(trait = trait, n_gwas_proteins = length(unique(sig_snp_genes)))
		return(out)
	})
	write.table(gwas_p_dat, file = n_gwas_proteins_file, 
			   row.names = F, col.names = T, quote = F, sep = "\t")
	rm(snp_genes)
}

gwas_comb_dat <- left_join(gwas_g_dat, gwas_p_dat)

generate_res_filenames <- function(df)
{
	pre_nam <- with(df, paste0("results/", trait, "/", trait, "_", sampling_method, "_sampling_", pd, "_"))
	out <- df %>%
		mutate(pathway_overlap_OR_filename = paste0(pre_nam, "pathway_enrichment_tests.txt"), 
			   enrichment_corr_filename = paste0(pre_nam, "pathway_enrichment_score_correlations.txt"), 
			   overlap_filename = paste0(pre_nam, "pathway_and_gene_overlap.txt"), 
			   gene_overlap_OR_filename = paste0(pre_nam, "gene_enrichment_tests.txt"))
	return(out)
}				




all_filenames <- list.files("results/glucose")

temp <- expand.grid(trait = traits,
 					pd = pathway_databases,
 					sampling_method = sampling_methods)

filenames_df <- generate_res_filenames(temp) %>% as_tibble


# to write out: 

# trait n-ewas-genes n-gwas-genes n-pathway-overlap observed-enrichment-correlation expected-enrichment-correlation p-diff

### get info from pathway_enrichment_score_correlations.txt

# gene tables:
# trait n-ewas-genes n-gwas-genes n-gene-overlap observed-overlap expected-overlap p-diff


# gene overlap results

z_test <- function(x, dist) {
	mu <- mean(dist)
	z <- (x - mu) / sd(dist)
	p <- 2 * pnorm(-abs(z))
	return(p)
}

# two tables needed one with all genes other with pcgs
filenames_df
gene_overlap_files <- filenames_df %>% 
	dplyr::filter(pd %in% c("go", "reactome")) %>%
	mutate(genes_or_proteins = case_when(pd == "go" ~ "genes", 
										 pd == "reactome" ~ "proteins"))

x=1
gene_overlap_tabs <- lapply(seq_along(unique(gene_overlap_files$pd)), function(x) {
	pathway_db <- unique(gene_overlap_files$pd)[x]
	temp <- gene_overlap_files %>%
		dplyr::filter(pd == pathway_db) %>%
		dplyr::filter(sampling_method == "random")
	out_tab <- map_dfr(1:nrow(temp), function(row) {
		print(row)
		df <- temp[row, ]
		if (df$pd %in% c("kegg", "go")) {
			n_gwas_val <- "n_gwas_genes"
		} else {
			n_gwas_val <- "n_gwas_proteins"
		}
		gene_overlap_out <- read_tsv(df$gene_overlap_OR_filename) %>%
				mutate(trait = df$trait, pd = "gene_only", sm = df$sampling_method) %>%
				dplyr::select(trait, pd, sm, everything())
		overlap_out <- read_tsv(df$overlap_filename) %>%
			mutate(trait = df$trait) %>%
			dplyr::select(trait, everything())
		comb_res <- gene_overlap_out %>%
			left_join(overlap_out) %>%
			left_join(gwas_comb_dat)
		# actual res
		a_res <- comb_res[comb_res$perm == 1, ]
		# permutation res
		p_res <- comb_res[comb_res$perm != 1, ]
		pval <- z_test(a_res$OR, p_res$OR)

		out <- a_res %>%
			dplyr::select(trait, n_ewas_genes, n_gwas_genes = one_of(n_gwas_val),
			 			  gene_overlap, OR) %>%
			rename(observed_OR = OR) %>%
			mutate(expected_gene_overlap = round(median(p_res$gene_overlap)), 
				   expected_OR = median(p_res$OR), 
				   p_diff = pval)
		return(out)
	})
	return(out_tab)
})
names(gene_overlap_tabs) <- unique(gene_overlap_files$pd)

## fix p-diff if expected OR is 0...
gene_overlap_tabs <- lapply(gene_overlap_tabs, function(x) {
	x %>% 
		mutate(p_diff = ifelse(expected_OR == 0 & observed_OR == 0, 1, p_diff))
})

pathway_db <- unique(filenames_df$pd)[1]
pathway_overlap_tabs <- lapply(unique(filenames_df$pd), function(pathway_db) {
	temp <- filenames_df %>%
		dplyr::filter(pd == pathway_db, 
				      sampling_method == "random")
	row=2
	out_tab <- map_dfr(1:nrow(temp), function(row) {
		df <- temp[row, ]
		if (df$pd %in% c("kegg", "go")) {
			n_gwas_val <- "n_gwas_genes"
		} else {
			n_gwas_val <- "n_gwas_proteins"
		}
		pathway_overlap_out <- read_tsv(df$pathway_overlap_OR_filename) %>%
				mutate(trait = df$trait, pd = pathway_db, sm = df$sampling_method) %>%
				dplyr::select(trait, pd, sm, everything())
		overlap_out <- read_tsv(df$overlap_filename) %>%
			mutate(trait = df$trait) %>%
			dplyr::select(trait, everything())
		comb_res <- pathway_overlap_out %>%
			left_join(overlap_out) %>%
			left_join(gwas_comb_dat)
		# actual res
		a_res <- comb_res[comb_res$perm == 1, ]
		# permutation res
		p_res <- comb_res %>%
			dplyr::filter(perm != 1, 
						  OR != Inf)
		if (nrow(p_res) < 500) {
			p_res$OR <- Inf
		}
		pval <- z_test(a_res$OR, p_res$OR)

		out <- a_res %>%
			dplyr::select(trait, n_ewas_genes, n_gwas_genes = one_of(n_gwas_val),
			 			  geneset_overlap = pathway_overlap, OR) %>%
			rename(observed_OR = OR) %>%
			mutate(expected_OR = mean(p_res$OR), 
				   p_diff = pval, n_perms = nrow(comb_res))
		return(out)
	})
	return(out_tab)
})
names(pathway_overlap_tabs) <- unique(filenames_df$pd)

# check permutations! So enrichment tests worked -- alll goooood!
n_perms <- map(pathway_overlap_tabs, function(x) dplyr::select(x, trait, n_perms))

pathway_db <- unique(filenames_df$pd)[1]
pathway_enrich_tabs <- lapply(unique(filenames_df$pd), function(pathway_db) {
	temp <- filenames_df %>%
		dplyr::filter(pd == pathway_db, 
				      sampling_method == "random")
	row=1
	out_tab <- map_dfr(1:nrow(temp), function(row) {
		df <- temp[row, ]
		if (df$pd %in% c("kegg", "go")) {
			n_gwas_val <- "n_gwas_genes"
		} else {
			n_gwas_val <- "n_gwas_proteins"
		}
		enrich_out <- read_tsv(df$enrichment_corr_filename) %>%
				mutate(trait = df$trait, pd = pathway_db, sm = df$sampling_method) %>%
				dplyr::select(trait, pd, sm, everything())
		overlap_out <- read_tsv(df$overlap_filename) %>%
			mutate(trait = df$trait) %>%
			dplyr::select(trait, everything())
		if (nrow(enrich_out) == 1) {
			overlap_out <- overlap_out[overlap_out$perm == 1, ]
			comb_res <- enrich_out %>%
				left_join(overlap_out) %>%
				left_join(gwas_comb_dat)
			out <- comb_res %>%
				mutate(expected_cor = 0, expected_geneset_overlap = 0) %>%
				dplyr::select(trait, n_ewas_genes, n_gwas_genes = one_of(n_gwas_val), 
							  observed_cor = cor_est, expected_cor, p_diff = cor_p)
			return(out)
		}
		comb_res <- enrich_out %>%
			left_join(overlap_out) %>%
			left_join(gwas_comb_dat)
		# actual res
		a_res <- comb_res[comb_res$perm == 1, ]
		# permutation res
		p_res <- comb_res[comb_res$perm != 1, ]
		pval <- z_test(a_res$cor, p_res$cor)

		out <- a_res %>%
			dplyr::select(trait, n_ewas_genes, n_gwas_genes = one_of(n_gwas_val), 
						  geneset_overlap = pathway_overlap, cor) %>%
			rename(observed_cor = cor) %>%
			mutate(expected_cor = mean(p_res$cor), p_diff = pval)
		return(out)
	})
	return(out_tab)
})
names(pathway_enrich_tabs) <- unique(filenames_df$pd)

save(gene_overlap_tabs, file = file.path(home_dir, "report/report_data/gene_overlap_results.RData"))
save(pathway_overlap_tabs, file = file.path(home_dir, "report/report_data/pathway_overlap_results.RData"))
save(pathway_enrich_tabs, file = file.path(home_dir, "report/report_data/pathway_enrichment_results.RData"))

# save(enrichment_res, file = "report/report_data/summarised_enrichment_results.RData")

