# ------------------------------------------------------
# Visualising and assessing correlation between enrichment scores
# ------------------------------------------------------

pkgs <- c("tidyverse", "gplots", "ggrepel", "gridExtra", "usefunc")
lapply(pkgs, require, character.only = T)

### read in all the data
opengwas_snps <- read_tsv("data/derived/ieugwasdb_top_hits.txt")
snp_genes <- read_tsv("data/derived/genetic/snps_mapped_to_genes.txt")

snp_genes <- snp_genes %>%
	dplyr::filter(name %in% opengwas_snps$SNP) %>%
	dplyr::select(-chromosome, -position, -hgnc_symbol, -distance_to_gene)

opengwas_snps <- opengwas_snps %>%
	dplyr::filter(SNP != "none") %>%
	left_join(snp_genes, by = c("SNP" = "name"))

opengwas_studies <- read_tsv("data/derived/ieugwasdb_studies.txt")

# duplicated gwas correlation res
duplicated_cor_res <- read_tsv("data/derived/duplicated_gwas_enrichment_score_correlations.txt")

# all correlations
correlations <- new_load("data/derived/all_enrichment_score_correlations.RData")
cor_test_res <- read_tsv("results/go_all_pathway_enrichment_correlations.txt")

# ewas traits
ewas_traits <- readLines("data/traits.txt")

# ------------------------------------------------------
# Visualise n-snps across traits
# ------------------------------------------------------

gene_counts <- opengwas_snps %>%
	dplyr::filter(SNP != "none") %>%
	dplyr::select(-SNP) %>%
	distinct() %>%
	group_by(id) %>%
	summarise(n_genes = n())

gene_counts %>%
	dplyr::filter(n_genes > 500)

n_gene_hist <- ggplot(gene_counts, aes(x = n_genes)) +
	geom_histogram(fill = "blue", colour = "black") +
	theme_bw()

ggsave("results/plots/n_genes_histogram.pdf", plot = n_gene_hist)

# gwas with fewer than 10 genes
sum(gene_counts$n_genes < 10) / nrow(gene_counts)
# 65% have fewer than 10 genes identified......

# ------------------------------------------------------
# Check correlations of GWAS of same traits
# ------------------------------------------------------

id_and_n <- dup_studies %>%
	dplyr::select(id, sample_size)

# checking snp count differences
duplicated_cor_res <- duplicated_cor_res %>%
	left_join(gene_counts, by = c("id1" = "id")) %>%
	rename(id1_ngene = n_genes) %>%
	left_join(gene_counts, by = c("id2" = "id")) %>%
	rename(id2_ngene = n_genes) %>%
	arrange(trait)

duplicated_cor_res %>%
	dplyr::filter(trait == "Neuroticism")

# would expect low concordance even when little difference in genes IF
# the total number of genes is low
# could also be a function of the trait. Would expect clinically measured
# traits to be the same, but mental health traits could be measured wildly
# differently across studies...

plot_res <- duplicated_cor_res %>%
	mutate(abs_gene_diff = abs(id1_ngene - id2_ngene), 
		   prop_gene_diff = case_when(id1_ngene < id2_ngene ~ id2_ngene / id1_ngene, 
		   							 id2_ngene < id1_ngene ~ id1_ngene / id2_ngene), 
		   total_ngene = id1_ngene + id2_ngene) 

# remove outlier res! 

plot_res %>%
	dplyr::filter(total_ngene > 1000) %>%
	pull(trait)

# separate out traits that are present 10 or more times
trait_n <- plot_res %>%
	group_by(trait) %>%
	summarise(count = n()) %>%
	arrange(desc(count))

trait_n <- trait_n %>%
	mutate(trait_count = 1:nrow(.)) %>%
	mutate(t_group = case_when(count < 10 ~ as.double(100), 
							   count > 9 ~ as.double(trait_count)))

# plot for each trait
x=trait_n$t_group[1]
p_list <- lapply(unique(trait_n$t_group), function(x) {
	traits <- trait_n %>%
		dplyr::filter(t_group == x) %>%
		pull(trait)
	p_res <- plot_res %>%
		dplyr::filter(trait %in% traits)
	p <- ggplot(p_res, aes(x = abs_gene_diff, y = rho, colour=total_ngene)) +
		geom_point() +
		theme_bw() + 
		ylim(0, 1)
	if (length(traits) > 1) {
		p <- p + 
			labs(title = "Other traits")
	} else {
		p <- p + 
			labs(title = traits)
	}
	return(p)
})
glist <- lapply(p_list, ggplotGrob)
pdf("results/plots/gene_diff_rho_plots.pdf") 
marrangeGrob(glist, nrow=1, ncol=1, top = NULL) # This takes forever for some reason...
dev.off()

# NOW FOR PROPORTIONAL DIFFERENCES! 
p_list <- lapply(unique(trait_n$t_group), function(x) {
	traits <- trait_n %>%
		dplyr::filter(t_group == x) %>%
		pull(trait)
	p_res <- plot_res %>%
		dplyr::filter(trait %in% traits) %>%
		dplyr::filter(trait != "Platelet count") %>%
		dplyr::filter(trait != "Rheumatoid arthritis")
	p <- ggplot(p_res, aes(x = prop_gene_diff, y = rho, colour=total_ngene)) +
		geom_point() +
		theme_bw() + 
		ylim(0, 1)
	if (length(traits) > 1) {
		p <- p + 
			geom_label_repel(aes(label = trait))
	} else {
		p <- p + 
			labs(title = traits)
	}
	return(p)
})
glist <- lapply(p_list, ggplotGrob)
pdf("results/plots/prop_gene_diff_rho_plots.pdf") 
marrangeGrob(glist, nrow=1, ncol=1, top = NULL) # This takes forever for some reason...
dev.off()

# ------------------------------------------------------
# Check histogram of rhos to see if z-test valid
# ------------------------------------------------------

rho_hist <- ggplot(cor_test_res, aes(x = rho)) +
	geom_histogram(colour = "black", fill = "blue") + 
	theme_bw()

ggsave("results/plots/rho_histogram.pdf", plot = rho_hist)

log_cor_test_res <- cor_test_res %>%
	mutate(log_rho = log(abs(rho)))

log_rho_hist <- ggplot(log_cor_test_res, aes(x = log_rho)) +
	geom_histogram(colour = "black", fill = "blue") + 
	theme_bw()

ggsave("results/plots/log_rho_histogram.pdf", plot = log_rho_hist)

summary(cor_test_res$rho)

# if big diff between mean and median then do a kruskal-wallis, otherwise
# z-test (in next section)

# ------------------------------------------------------
# EWAS vs. all GWAS
# ------------------------------------------------------

# calculate pval threshold
ncor <- nrow(cor_test_res)
bon_threshold <- 0.05 / ncor
# ewas-only p value threshold
# ewas_ncor <- length(ewas_traits) * nrow(correlations) # roughly
# ewas_threshold <- 0.05 / ewas_ncor
ewas_threshold <- bon_threshold

# separate out ewas dat
ewas_only_dat <- cor_test_res %>%
	dplyr::filter(trait1 %in% ewas_traits & trait2 %in% ewas_traits)

ewas_gwas_dat <- cor_test_res %>%
	dplyr::filter(trait1 %in% ewas_traits | trait2 %in% ewas_traits) %>%
	dplyr::filter(!(trait1 %in% ewas_traits & trait2 %in% ewas_traits))

gwas_only_dat <- cor_test_res %>%
	dplyr::filter(!(trait1 %in% ewas_traits) & !(trait2 %in% ewas_traits))

# compare rho to the average 
z_test <- function(x, dist) {
	mu <- mean(dist)
	z <- (x - mu) / sd(dist)
	p <- 2 * pnorm(-abs(z))
	return(p)
}

ewas_only_dat <- ewas_only_dat %>%
	mutate(p_diff = z_test(rho, cor_test_res$rho)) %>%
	mutate(p_diff_adj = p.adjust(p_diff, method = "BH", n = ncor))

ewas_gwas_dat <- ewas_gwas_dat %>%
	mutate(p_diff = z_test(rho, cor_test_res$rho)) %>%
	mutate(p_diff_adj = p.adjust(p_diff, method = "BH", n = ncor))

sum(ewas_gwas_dat$p_diff_adj < 0.05)
ewas_gwas_dat %>%
	dplyr::filter(p_diff < 0.05)

gwas_only_dat <- gwas_only_dat %>%
	mutate(p_diff = z_test(rho, cor_test_res$rho)) %>%
	mutate(p_diff_adj = p.adjust(p_diff, method = "BH", n = ncor))

gwas_out_dat <- gwas_only_dat %>%
	dplyr::filter(p_diff_adj < 0.05) %>%
	arrange(desc(rho))

mean(ewas_gwas_dat$rho)

# ------------------------------------------------
# Things to output
# ------------------------------------------------

gwas_studies <- opengwas_studies %>%
	dplyr::filter(tolower(id) %in% tolower(colnames(correlations)))

grep("lung", gwas_studies$trait, ignore.case = T, value = T)
lc_id <- gwas_studies %>%
	dplyr::filter(trait == "Lung cancer") %>%
	pull(id)

lc_sm_cor <- cor_test_res %>%
	dplyr::filter(trait1 == "current_versus_never_smoking", trait2 == lc_id)

out_list <- list(cor_range = range(cor_test_res$rho), 
				 cor_mean = mean(cor_test_res$rho),
				 ewas_gwas_mean = mean(ewas_only_dat$rho), 
				 gwas_gwas_mean = mean(gwas_only_dat$rho), 
				 gwas_cor_tab = gwas_out_dat, 
				 gwas_studies = gwas_studies,
				 lc_sm_cor = lc_sm_cor)

save(out_list, file="report/report_data/ewas_gwas_comp_stats.RData")

# range of correlations
# mean of all correlations
# mean of ewas-gwas correlations
# mean of gwas-gwas correlations
# table of gwas-gwas correlations greater than expected by chance
# gwas studies table with only traits that had SNPs at P<5e-8!

# ------------------------------------------------
# clustering stuuuuff
# ------------------------------------------------
D <-  as.dist(1 - abs(correlations))
PRhoTree <- hclust(D, method = "complete")
k <- cutree(PRhoTree, h = 0.2) ## distances are in 1 - Pearson Rho distances
table(k) ## 42 Clusters

# making the clusters
hc_clusters <- tibble(label = names(k), clust = k)
clusters <- count(hc_clusters, clust) %>%
	dplyr::filter(n > 1) %>%
	pull(clust)

out_clusters <- hc_clusters %>%
	dplyr::filter(clust %in% clusters) %>%
	arrange(clust)

clustnam <- paste0("results/go_clusters.txt")
write.table(out_clusters, file = clustnam, 
			col.names = T, row.names = F, quote = F, sep = "\t")

clustnam <- paste0("report/report_data/go_clusters.txt")
write.table(out_clusters, file = clustnam, 
			col.names = T, row.names = F, quote = F, sep = "\t")
