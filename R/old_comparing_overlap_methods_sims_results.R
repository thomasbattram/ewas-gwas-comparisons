# -------------------------------------------------------
# Examining the simulation results
# -------------------------------------------------------

pkgs <- c("tidyverse", "gridExtra", "pROC", "RColorBrewer")
lapply(pkgs, require, character.only = T)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

devtools::load_all("~/repos/usefunc/")

# -------------------------------------------------------
# Read in data and bind it together
# -------------------------------------------------------
bind_data <- function(file_dir, file_prefix) 
{
	# extracts data from simulation results

	files <- list.files(file_dir)
	files <- grep(file_prefix, files, value = TRUE)
	if (length(files) == 0) stop("Files are not in the directory provided")
	# extract data in loop and bind it together
	file_paths <- file.path(file_dir, files)
	out <- map_dfr(file_paths, read_tsv)
	return(out)
}

temp_dir <- "data/derived/sim_temp"
file_prefix <- "method_sim_"
res_file <- "results/sim/method_sim_res.txt"
if (file.exists(res_file)) {
	sim_res <- read_tsv(res_file)
	sim_n <- length(unique(sim_res$sim))
	if (sim_n != 1000) {
		sim_res <- bind_data(temp_dir, file_prefix)
		write.table(sim_res, file = res_file, 
					col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	}
} else {
	sim_res <- bind_data(temp_dir, file_prefix)
	write.table(sim_res, file = res_file, 
				col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
}

#
# old_res <- read_tsv("results/sim/old_method_sim_res.txt") %>%
# 	dplyr::select(-or_p, -or_g, -pathway_overlap, -gene_overlap)

# sim_res <- sim_res %>%
# 	dplyr::filter(rho_p == 1) %>%
# 	dplyr::select(-rho_p) %>%
# 	left_join(old_res)

# write.table(sim_res, file = res_file, 
# 			col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# -------------------------------------------------------
# Check results!
# -------------------------------------------------------

roc_res1 <- with(sim_res, roc(percent_ewas_causal ~ or_p))

roc_res2 <- with(sim_res, roc(percent_ewas_causal ~ rho_p))

roc_res3 <- with(sim_res, roc(percent_ewas_causal ~ or_g))

# ewas power and ngenes and pathway_db...

### overall predictive power in each pathway across all vars
unique_pdb <- unique(sim_res$pathway_db)
variables <- c("or_g", "or_p", "rho_p")
colours <- brewer.pal(n = length(variables), name = "Set1")

overall_res <- lapply(unique_pdb, function(pdb) {
	sim_res_temp <- sim_res %>%
		dplyr::filter(pathway_db == pdb)
	roc_list <- lapply(variables, function(x) {
		roc_out <- roc(sim_res_temp[["percent_ewas_causal"]] ~ sim_res_temp[[x]])
		return(roc_out)
	})
	names(roc_list) <- variables
	p_roc <- pROC::ggroc(roc_list) + 
		geom_abline(intercept = 1, slope = 1, colour = "black", alpha = 0.6) +
		labs(title = pdb) +
		scale_colour_manual(values = colours) +
		theme_bw()
	return(p_roc)
})
glist <- lapply(overall_res, ggplotGrob)
pdf("results/plots/methods_simulations_overall_res.pdf")
marrangeGrob(glist, nrow=1, ncol=1, top = NULL)
dev.off()

### res split by ewas power and ngenes -- just extract AUCs here!

extract_auc_ci <- function(roc_out) {
	auc <- t(as.data.frame(ci.auc(roc_out)))
	rownames(auc) <- NULL
	colnames(auc) <- c("lower_ci", "estimate", "upper_ci")
	auc <- as.data.frame(auc)
	return(auc)
}

get_auc <- function(var1, var2) {
	roc_out <- roc(var1 ~ var2)
	auc <- extract_auc_ci(roc_out)
	return(auc)
}

auc_res <- sim_res %>% 
	pivot_longer(
		cols = c("or_g", "or_p", "rho_p"),
		names_to = c("method"), 
		values_to = "value"
	) %>%
	group_by(pathway_db, ewas_power, n_genes, method) %>%
	summarise(auc_lower_ci = get_auc(percent_ewas_causal, value)$lower_ci, 
			  auc_estimate = get_auc(percent_ewas_causal, value)$estimate,
			  auc_upper_ci = get_auc(percent_ewas_causal, value)$upper_ci)

auc_plot <- ggplot(auc_res, aes(x = as.factor(ewas_power), y = auc_estimate, colour = method)) +
				geom_point(position = position_dodge(width = 0.9)) +
				geom_linerange(aes(ymin = auc_lower_ci, ymax = auc_upper_ci), position = position_dodge(width = 0.9)) +
				labs(x = "ewas power", y = "AUC") +
				facet_grid(n_genes ~ pathway_db) +
				theme_bw()

ggsave('results/plots/methods_test_auc_plot.pdf', plot = auc_plot)

# checking association between actual gene overlap and ORs/Rhos

gene_o_res <- sim_res %>% 
	mutate(or_g = log(or_g), or_p = log(or_p)) %>%
	pivot_longer(
		cols = c("or_g", "or_p", "rho_p"),
		names_to = c("method"), 
		values_to = "value"
	) %>%
	group_by(gene_overlap, gwas_power, method, pathway_db, n_genes, ewas_power) %>%
	summarise(med_val = median(value))


variables <- c("or_g", "or_p", "rho_p")
gene_o_plots <- lapply(variables, function(x) {
	df <- gene_o_res %>%
		dplyr::filter(method == x)
	gene_o_plot <- 	ggplot(df, aes(x = gene_overlap, y = med_val)) +
						geom_point() + 
						# geom_line() + 
						facet_grid(n_genes ~ pathway_db) + 
						labs(title = x)
						theme_bw()
	return(gene_o_plot)
})
glist <- lapply(gene_o_plots, ggplotGrob)
pdf("results/plots/methods_simulations_gene_overlap.pdf")
marrangeGrob(glist, nrow=1, ncol=1, top = NULL)
dev.off()


gene_o_plot <- ggplot(gene_o_res, aes(x = gene_overlap, y = value, colour = method)) +
	geom_point() + 
	geom_line() + 
	facet_grid(n_genes ~ pathway_db) + 
	theme_bw()

ggsave("test.pdf", gene_o_plot)















