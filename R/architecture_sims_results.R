# -------------------------------------------------------
# Examining the simulation results
# -------------------------------------------------------

## takes sim results from architecture_sims.R, binds them together,
## assesses whether there are any scenarios (overlap between causal and
## associated genes along with number of genes to discover)

pkgs <- c("tidyverse", "gridExtra", "RColorBrewer", "pROC", "usefunc")
lapply(pkgs, require, character.only = T)

library(cowplot) # to arrange plots nicely
library(egg)

# setwd("SCRATCH_SPACE")
home_dir <- "~/projects/epi_gen_comp" # CHANGE ME WHEN NEEDED

# -------------------------------------------------------
# bind results together
# -------------------------------------------------------

sim_res_path <- "data/derived/sim_temp/architecture_sims/"

bind_results <- function(sim_file_path, res_dir, sim_nam) 
{
	### read in and bind results from simulations
	all_files <- list.files(sim_file_path)
	# sim_nam <- paste0("sim", sim_n, "_", data_used)
	sim_file_names <- paste0(sim_nam, "_", 1:100, ".txt")
	sim_files <- all_files[all_files %in% sim_file_names]
	message("Successful iterations: ", length(sim_files))
	out <- map_dfr(sim_files, function(sf) {
		suppressMessages(dat <- read_tsv(file.path(sim_file_path, sf)))
		message("returning data")
		return(dat)
	})
	out_nam <- file.path(res_dir, paste0(sim_nam, ".txt"))
	message("Saving results to ", out_nam)
	write.table(out, out_nam, 
				col.names = T, row.names = F, quote = F, sep = "\t")
	return(NULL)
}

res_dir <- "results/sim"
bind_results(sim_file_path = sim_res_path, res_dir = res_dir, 
			 sim_nam = "architecture_sims_go")

# -------------------------------------------------------
# read in the results
# -------------------------------------------------------
traits <- readLines("data/traits.txt")
res_dir <- "results/sim"
sim_res <- map_dfr(traits, function(trait) {
	res_nam <- paste0("architecture_sims_go_", trait, ".txt")
	file_path <- file.path(res_dir, res_nam)
	if (!file.exists(file_path)) {
		warning(file_path, " does not exist")
		return(NULL)
	}
	read_tsv(file_path)
})

## swap out glucose + insulin for fasting gluc + ins
sim_res$trait <- gsub("insulin", "fasting_insulin", sim_res$trait)
sim_res$trait <- gsub("glucose", "fasting_glucose", sim_res$trait)

# -------------------------------------------------------
# get empirical data for traits
# -------------------------------------------------------

## for tests: traits <- traits[traits != "c-reactive_protein"]
pathway_databases <- c("gene_ontology")
sampling_methods <- c("random")
temp <- expand.grid(trait = traits,
 					pd = pathway_databases,
 					sampling_method = sampling_methods)
row <- 1
trait_res <- map_dfr(1:nrow(temp), function(row) {
	df <- temp[row, ]
	df$pd <- ifelse(df$pd == "gene_ontology", "go", "kegg")
	pre_nam <- paste0("results/", df$trait, "/",
					  df$trait, "_", df$sampling_method, "_sampling_",  df$pd, "_")
	enrich_nam <- paste0(pre_nam, "pathway_enrichment_tests.txt")
	enrich_out <- read_tsv(enrich_nam) %>%
		mutate(trait = df$trait, or_p = OR) %>%
		dplyr::select(trait, perm, or_p)
	overlap_nam <- paste0(pre_nam, "pathway_and_gene_overlap.txt")
	overlap_out <- read_tsv(overlap_nam) %>%
		mutate(trait = df$trait)
	gene_enrich_nam <- paste0(pre_nam, "gene_enrichment_tests.txt")
	gene_enrich <- read_tsv(gene_enrich_nam) %>%
		mutate(trait = df$trait, or_g = OR) %>%
		dplyr::select(trait, perm, or_g)
	rho_enrich_nam <- paste0(pre_nam, "pathway_enrichment_score_correlations.txt")
	rho_enrich_out <- read_tsv(rho_enrich_nam) %>%
		mutate(trait = df$trait, rho_p = cor) %>%
		dplyr::select(trait, perm, rho_p)
	comb_res <- enrich_out %>%
		left_join(gene_enrich) %>%
		left_join(rho_enrich_out) %>%
		left_join(overlap_out) %>%
		dplyr::filter(perm == 1)
	return(comb_res)
})

gwas_genes_file <- "data/derived/gwas_genes.txt"
gwas_g_dat <- read_tsv(gwas_genes_file)
gwas_ng_dat <- gwas_g_dat %>%
	group_by(trait) %>%
	summarise(n_gwas_genes = n())

trait_res <- trait_res %>%
	left_join(gwas_ng_dat)

# -------------------------------------------------------
# Assess whether empirical results are different to the simulations
# -------------------------------------------------------

z_test <- function(x, dist) {
	mu <- mean(dist, na.rm=T)
	z <- (x - mu) / sd(dist, na.rm=T)
	p <- 2 * pnorm(-abs(z))
	return(p)
}

# could do a z-test for difference at each thingy
trait_res2 <- trait_res %>%
	dplyr::select(trait, true_cor = rho_p, true_or_g = or_g, 
				  n_gwas_genes, n_ewas_genes)

z_test_res <- sim_res %>%
	left_join(trait_res2) %>% 
	mutate(N_KTG = n_gwas_genes + n_ewas_genes, 
		   R_KTG = paste0(n_gene_proportions, ":1"), 
		   percentage_discovered = 1 / (n_gene_proportions*2) * 100) %>%
	group_by(trait, N_KTG, R_KTG, percentage_discovered, ca_con_overlap) %>%
	summarise(mean_sim_cor = mean(rho_p, na.rm=T), 
			  mean_sim_or_g = mean(or_g),
			  p_diff_cor = z_test(unique(true_cor), rho_p), 
			  p_diff_or_g = z_test(unique(true_or_g), or_g)) %>%
	mutate(p_diff_cor_fdr = p.adjust(p_diff_cor, method = "bonferroni"), 
		   p_diff_or_g_fdr = p.adjust(p_diff_or_g, method = "bonferroni"))

# check p diffs
gene_vars <- grep("or_g", c(colnames(z_test_res), colnames(trait_res2)), value = T)
cor_vars <- grep("_cor", c(colnames(z_test_res), colnames(trait_res2)), value = T)

# correlation
cor_z_test_res <- z_test_res %>%
	dplyr::filter(p_diff_cor_fdr < 0.05) %>%
	left_join(trait_res2) %>%
	dplyr::select(-one_of(gene_vars), -p_diff_cor_fdr) %>%
	dplyr::select(trait, N_KTG, R_KTG, percentage_discovered, ca_assoc_overlap = ca_con_overlap,
				  observed_cor = true_cor, expected_cor = mean_sim_cor, 
				  p_diff = p_diff_cor)

write.table(cor_z_test_res, file = file.path(home_dir, "report/report_data/arch_sims_ztest_sig_res.tsv"), 
			col.names = T, row.names = F, sep = "\t")

# -------------------------------------------------------
# Plot the geneset enrichment score simulation and empirical results
# -------------------------------------------------------

### FUNCTIONS FOR EVERYTHING BABY!
extract_sig_z_test <- function(z_dat, trait_nam)
{
	z_dat %>%
		ungroup() %>%
		dplyr::filter(trait == trait_nam) %>%
		mutate(n_gene_proportions = as.numeric(gsub(":1", "", R_KTG))) %>%
		mutate(to_colour = TRUE) %>%
		dplyr::select(trait, ca_con_overlap = ca_assoc_overlap, 
					  n_gene_proportions, to_colour)
}

extract_empirical_rho <- function(trait_dat, trait_nam)
{
	trait_dat %>%
		dplyr::filter(trait == trait_nam) %>%
		dplyr::rename(OR = one_of("rho_p"))		
}

sort_sim_res_for_plot <- function(sim_dat, trait_nam, z_dat) 
{
	out <- sim_dat %>%
		dplyr::filter(trait == trait_nam) %>%
		dplyr::rename(OR = one_of("rho_p")) %>%
		left_join(z_dat) %>%
		mutate(to_colour = ifelse(is.na(to_colour), FALSE, to_colour))
	# out$n_gene_proportions2 <- factor(paste0(out$n_gene_proportions, ":1"), levels = c("1:1", "2:1", "3:1", "5:1", "10:1", "20:1"))
	out$n_gene_proportions2 <- factor(paste0("1/", out$n_gene_proportions * 2), levels = c("1/2", "1/4", "1/6", "1/10", "1/20", "1/40"))
	return(out)
}

plot_sim_dat <- function(sim_dat, colours, trait_dat, trait_nam, two_traits = FALSE)
{
	# make the number of genes text for plot
	# neg <- paste0("N[KEG]: ", trait_dat$n_ewas_genes)
	# ngg <- paste0("N[KGG]: ", trait_dat$n_gwas_genes)
	# ntg <- paste0("N[KTG]: ", trait_dat$n_ewas_genes + trait_dat$n_gwas_genes)

	text_size <- ifelse(two_traits, 5, 3)

	p <- ggplot(sim_dat) +
		geom_boxplot(aes(x = as.factor(ca_con_overlap), y = OR,
						 fill = n_gene_proportions2, colour = to_colour)) +
		# scale_colour_discrete(name = bquote("Ratio of causal and associated genes relative to " ~ N[KTG])) + 
		scale_fill_manual(name = "Fraction of genes GWAS and EWAS have discovered        ", # have to put loads of spaces afterwards because the plotting function is dumb - try it without
						    values = cols_to_use) +
		scale_colour_manual(name = "",
						  values = c("TRUE" = "grey", "FALSE" = "black")) +
		guides(colour = "none") + 
		labs(x = "Proportion of simulated associated genes that are causal",
			 y = "Geneset enrichment score correlation") +
		ylim(0, 1) + 
		# annotate("text", x = 1, y = 0.8, label = neg, parse = TRUE, size = text_size) + 
		# annotate("text", x = 2.25, y = 0.8, label = ngg, parse = TRUE, size = text_size) + 
		# annotate("text", x = 3.5, y = 0.8, label = ntg, parse = TRUE, size = text_size) + 
		geom_hline(yintercept = trait_dat$OR, colour = "red", linetype = "dashed") +
		labs(title = str_to_title(gsub("_", " ", trait_nam))) + 
		theme_bw() + 
		theme(plot.title = element_text(margin = margin(b = -15))) 
	
	if (!two_traits) {
		p <- p +
			theme(text = element_text(size = 7))
	}
	
	return(p)
}

#cols_to_use <- brewer.pal(n = 6, name = 'Dark2')
cols_to_use <- get_cb_palette()[1:6]
names(cols_to_use) <- c("1/2", "1/4", "1/6", "1/10", "1/20", "1/40")

## just crp and urate
crp_fvns <- c("c-reactive_protein", "former_versus_never_smoking")
crp_fvns_plots <- lapply(seq_along(crp_fvns), function(x) {
	tra <- crp_fvns[x]
	print(tra)	

	# get correlation of geneset enrichment scores from empirical analyses
	trait_dat <- extract_empirical_rho(trait_res, tra)

	# sort evidence of differences between empirical and simulated enrichment 
	# score correlations for use in plot
	z_test_dat <- extract_sig_z_test(cor_z_test_res, tra)

	# gather simulation data and bind with z-test data for plot
	sim_dat <- sort_sim_res_for_plot(sim_res, tra, z_test_dat)

	# make boxplots 
	p_out <- plot_sim_dat(sim_dat, cols_to_use, trait_dat, tra, two_traits = TRUE)

# ggsave("results/plots/test.pdf", plot = p_out)

	return(p_out)
})
names(crp_fvns_plots) <- crp_fvns

# ## get legend
# legend <- get_legend(
# 	crp_fvns_plots[[1]] + theme(legend.position="bottom")
# )

## set x and y axes
new_crp_fvns_plots <- lapply(seq_along(crp_fvns_plots), function(x) {
	if (x == length(crp_fvns_plots)) {
		out <- crp_fvns_plots[[x]]
	} else {
	out <- crp_fvns_plots[[x]] + 
		   		theme(axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.x = element_blank())
	}
	return(out)
})

outp <- ggpubr::ggarrange(plotlist = new_crp_fvns_plots, nrow = 2, ncol = 1, common.legend = TRUE, legend="bottom")

## output plot
# outp <- plot_grid(crp_fvns_plots[["c-reactive_protein"]] + theme(legend.position="none"),
# 				  crp_fvns_plots[["former_versus_never_smoking"]] + theme(legend.position="none"),
# 				  # label_x = 1, label_y = 1,
# 				  hjust = -0.3,
# 				  ncol = 1, nrow = 2)

# outp <- plot_grid(new_crp_fvns_plots[[1]], legend, ncol=1, nrow=2, rel_heights = c(1, 0.1))

ggsave("results/plots/architecture_sims_crp_fvns_only_correlation_of_pathway_enrichment_scores.png", 
	   plot = outp, bg = "transparent")


### now other traits!
x=8
pl <- lapply(seq_along(traits), function(x) {
	tra <- traits[x]
	print(tra)	

	# get correlation of geneset enrichment scores from empirical analyses
	trait_dat <- extract_empirical_rho(trait_res, tra)

	# sort evidence of differences between empirical and simulated enrichment 
	# score correlations for use in plot
	z_test_dat <- extract_sig_z_test(cor_z_test_res, tra)

	# gather simulation data and bind with z-test data for plot
	sim_dat <- sort_sim_res_for_plot(sim_res, tra, z_test_dat)

	# make boxplots 
	p_out <- plot_sim_dat(sim_dat, cols_to_use, trait_dat, tra)

# ggsave("results/plots/test.pdf", plot = p_out)

	return(p_out)
})


names(pl) <- traits

## Remove traits used in main plot
pl <- pl[!names(pl) %in% crp_fvns]

## Make sure to use a trait here that will have all 6 levels for the legend
legend <- get_legend(
	pl[["fasting_glucose"]] + theme(legend.position = "bottom")
)

pl <- lapply(pl, function(x) {x + theme(legend.position = "none")})

bottom_pl <- names(pl)[length(pl)]
bottom_row <- plot_grid(pl[[bottom_pl]], legend, nrow=1, rel_widths = c(1, 2))

outp <- plot_grid(plotlist = pl[names(pl) != bottom_pl], ncol = 3, nrow = 4)

# outp <- plot_grid(outp, legend, ncol=1, nrow=2, rel_heights = c(1, 0.1))
outp <- plot_grid(outp, bottom_row, ncol=1, nrow=2, rel_heights = c(4, 1))

ggsave("results/plots/architecture_sims_other_traits_corr.png", plot = outp, width = 12, height = 10, 
	   bg = "transparent")

# gene overlap OR
z_test_res %>%
	dplyr::filter(p_diff_or_g_fdr < 0.05) %>%
	left_join(trait_res2) %>%
	dplyr::select(-one_of(cor_vars), -p_diff_or_g_fdr)


# where ca_con_overlap = 0, it's saying that the only causal genes the 
# ewas is picking up are those that have overlapped with the gwas
# and if this was the case, what would we expect overlap to be 
