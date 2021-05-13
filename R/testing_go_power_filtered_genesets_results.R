# -------------------------------------------------------
# Examining the simulation results -- testing go term filtering
# -------------------------------------------------------

## This script aims to assess whether there is any difference in power
## when assessing the overlap of genes + genesets when limiting go term
## genesets to those with more than 1 and fewer than 5000 genes. It comes after
## the simulations are performed.

## pkgs
library(tidyverse) # tidy code + data
library(RColorBrewer) # colours for plots
library(pROC) # generating ROC curves and calculating AUC

devtools::load_all("~/repos/usefunc/") # misc useful functions

## directory
dir <- "~/main_project/epi_gen_comp/"
setwd(dir)


## structure:
# 1. bind new simulation results together and save them
# 2. read in old simulation results and bind to new ones
# 3. plot power as was done previously

# -------------------------------------------------------
# bind new simulation results together
# -------------------------------------------------------

sim_res_path <- "data/derived/sim_temp/geneset_filtering_tests"

bind_results <- function(sim_file_path, res_dir, sim_nam) 
{
	### read in and bind results from simulations
	all_files <- list.files(sim_file_path)
	sim_file_names <- paste0(sim_nam, "_", 1:100, ".txt")
	sim_files <- all_files[all_files %in% sim_file_names]
	message("Successful iterations: ", length(sim_files))
	out <- map_dfr(sim_files, function(sf) {
		suppressMessages(dat <- read_tsv(file.path(sim_file_path, sf)))
		if (!0 %in% dat$percent_ewas_causal) return(NULL)
		message("returning data")
		return(dat)
	})
	out_nam <- file.path(res_dir, paste0(sim_nam, ".txt"))
	message("Saving results to ", out_nam)
	write.table(out, out_nam, 
				col.names = T, row.names = F, quote = F, sep = "\t")
	return(NULL)
}

bind_results(sim_res_path, res_dir = "results/sim", sim_nam = "filter_sim")

# -------------------------------------------------------
# read in old simulation results and bind to new ones
# -------------------------------------------------------

old_res <- read_tsv("results/sim/sim2.txt", guess_max = 1e6)
new_res <- read_tsv("results/sim/filter_sim.txt", guess_max = 1e6)

## remove parameter values not present in new analyses
old_res <- old_res %>%
	dplyr::filter(pathway_type == "go", 
				  n_genes %in% c(500, 10000), 
				  gwas_power == 0.1, 
				  ewas_power == 0.2, 
				  percent_ewas_causal %in% c(0, 0.1) 
				  ) %>%
	mutate(pathway_type = paste0(pathway_type, "_non_filtered"))

## for tests
# new_res <- old_res

comb_res <- new_res %>%
	mutate(pathway_type = paste0(pathway_type, "_filtered")) %>%
	bind_rows(old_res)

# -------------------------------------------------------
# AUC and plotting functions
# -------------------------------------------------------

extract_auc_ci <- function(roc_out) 
{
	# extracts AUC and confidence intervals from ROC curve

	auc <- t(as.data.frame(ci.auc(roc_out)))
	rownames(auc) <- NULL
	colnames(auc) <- c("lower_ci", "estimate", "upper_ci")
	auc <- as.data.frame(auc)
	return(auc)
}

get_auc <- function(var1, var2) 
{
	# generates ROC curve and extracts AUC

	roc_out <- roc(var1 ~ var2)
	auc <- extract_auc_ci(roc_out)
	return(auc)
}

gen_auc_res <- function(sim_res) 
{
	### make ROC curves, get AUC
	parameters <- c("pathway_type", "ewas_power", "n_genes", "method")

	auc_res <- sim_res %>% 
		pivot_longer(
			cols = c("or_g", "or_p", "rho_p"),
			names_to = c("method"), 
			values_to = "value"
		) %>%
		group_by_at(parameters) %>%
		dplyr::summarise(auc_lower_ci = get_auc(percent_ewas_causal, value)$lower_ci, 
				  auc_estimate = get_auc(percent_ewas_causal, value)$estimate,
				  auc_upper_ci = get_auc(percent_ewas_causal, value)$upper_ci)
	return(auc_res)
}

plot_auc_res <- function(auc_res, percent_ewas_causal, title_nam) 
{
	### take auc res from gen_auc_res and plot them

	auc_res <- auc_res %>%
		mutate(method = case_when(method == "or_g" ~ "gene overlap", 
								  method == "or_p" ~ "geneset overlap", 
								  method == "rho_p" ~ "geneset enrichment correlation"))

	auc_plot <- ggplot(auc_res, aes(x = method, y = auc_estimate, 
								    colour = pathway_type, group = pathway_type)) +
				geom_point(position = position_dodge(width = 0.9)) +
				# geom_linerange(aes(ymin = auc_lower_ci, ymax = auc_upper_ci), position = "dodge") +
				geom_errorbar(aes(ymin = auc_lower_ci, ymax= auc_upper_ci), position = position_dodge(width = 0.9)) + 
				# geom_line() +
				# ylim(0.5, 1) +
				scale_y_continuous(breaks = seq(0.5, 1, 0.1), limits=c(0.47,1)) +
				labs(x = "Overlap method", y = "AUC", title = title_nam) + 
				theme_bw() + 
				# theme(axis.text.x = element_text(angle = 45)) +  
				facet_grid(n_genes ~ ., 
					   labeller = labeller(n_genes = n_genes_labs))		
	return(auc_plot)
}


# -------------------------------------------------------
# plot results 
# -------------------------------------------------------

## remember EWAS power == 0.2

auc_res <- gen_auc_res(comb_res)

n_genes_labs <- paste("n_genes:", unique(comb_res$n_genes))
names(n_genes_labs) <- as.character(unique(comb_res$n_genes))

auc_plot_title <- "The proportion of causal EWAS genes = 0.1"

auc_plot_out <- plot_auc_res(auc_res = auc_res, 
							 percent_ewas_causal = 0.1, 
							 title_nam = auc_plot_title)

ggsave("results/plots/go_sim1_filtered_terms_power.pdf", plot = auc_plot_out)

