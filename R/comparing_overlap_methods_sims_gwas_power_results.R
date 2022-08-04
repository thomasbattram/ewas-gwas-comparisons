# -------------------------------------------------------
# Examining the simulation results
# -------------------------------------------------------

## script takes simulation results from comparing_overlap_methods_sims.R
## binds them all together and examines the results 3 fold:
## 1. (main reason) to establish the best method and expected power
## 	  to detect gene/geneset overlap depending on scenario
## 2. assess whether there is any difference between mapping to all 
## 	  ensembl gene IDs or just protein coding genes when using
##	  GO or KEGG terms for geneset enrichment analyses

## NB. removed simulation 1 results as the "pathway down" method of
## extracting genes did not work. 

# also script is pretty messy... 

pkgs <- c("tidyverse", "gridExtra", "RColorBrewer", "pROC", "usefunc", "latex2exp")
lapply(pkgs, require, character.only = T)

# setwd(SCRATCH_SPACE)

sim_res_path <- "data/derived/sim_temp_gwas_power"

# -------------------------------------------------------
# bind results together
# -------------------------------------------------------

bind_results <- function(sim_n, sim_file_path, res_dir, sim_nam) 
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

sim_res_dir <- "results/sim_gwas_power"
sim_nam <- "sim2_go_genes"
bind_results(2, sim_res_path, res_dir = sim_res_dir, sim_nam = sim_nam)

## read in res 
sim_res_bound <- read_tsv(file.path(sim_res_dir, paste0(sim_nam, ".txt")))

# -------------------------------------------------------
# make functions to get auc of roc curves and plot them
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

gen_auc_res <- function(sim_res, sim_n) 
{
	### make ROC curves, get AUC
	if (sim_n == 1) {
		parameters <- c("pathway_type", "gwas_power", "prop_consequent_pathways", "method")
	} else if (sim_n == 2) {
		parameters <- c("pathway_type", "gwas_power", "n_genes", "method")
	} else {
		stop("sim_n should be 1 or 2")
	}

	auc_res <- sim_res %>% 
		pivot_longer(
			cols = c("or_g", "or_p", "rho_p"),
			names_to = c("method"), 
			values_to = "value"
		) %>%
		mutate(value = ifelse(value == Inf, 1e6, value)) %>%
		group_by_at(parameters) %>%
		dplyr::summarise(auc_lower_ci = get_auc(percent_ewas_causal, value)$lower_ci, 
				  auc_estimate = get_auc(percent_ewas_causal, value)$estimate,
				  auc_upper_ci = get_auc(percent_ewas_causal, value)$upper_ci)
	return(auc_res)
}

plot_auc_res <- function(auc_res, sim_n) 
{
	### take auc res from gen_auc_res and plot them

	auc_plot <- ggplot(auc_res, aes(x = as.factor(gwas_power), y = auc_estimate, 
								    colour = method, group = method)) +
				geom_point() +
				geom_linerange(aes(ymin = auc_lower_ci, ymax = auc_upper_ci)) +
				geom_line() +
				# ylim(0.5, 1) +
				scale_y_continuous(breaks = seq(0.5, 1, 0.1), limits=c(0.47,1)) +
				scale_colour_manual(name = "Method", 
									values = get_cb_palette()[1:3], 
									breaks = c("or_g", "or_p", "rho_p"), 
									labels = list(TeX(r'($OR_g$)'), TeX(r'($OR_p$)'), TeX(r'($\rho_p$)'))) +
				labs(x = "GWAS power", y = "AUC") + 
				theme_bw() # + 
				# theme(axis.text.x = element_text(angle = 45))
	if (sim_n == 1) {
		auc_plot <- auc_plot + 
			facet_grid(prop_consequent_pathways ~ pathway_type, 
					   labeller = labeller(prop_consequent_pathways = prop_consequent_pathways_labs))
	} else if (sim_n == 2) {
		auc_plot <- auc_plot + 
			facet_grid(n_genes ~ ., labeller = labeller(n_genes = n_genes_labs))
	} else {
		stop("sim_n should be 1 or 2")
	}
	return(auc_plot)
}

# -------------------------------------------------------
# Analyse the results for sim2 -- protein mapping
# -------------------------------------------------------

auc_res <- gen_auc_res(sim_res_bound, 2)

n_genes_labs <- paste("Assoc genes:", unique(auc_res$n_genes))
names(n_genes_labs) <- as.character(unique(auc_res$n_genes))

out_plot <- plot_auc_res(auc_res, 2)	
ggsave("results/plots/methods_test_gene_up_auc_plot_gwas_power.png", plot=out_plot)

auc_plot_list <- lapply(seq_along(auc_res_list), function(i) {
	pec <- names(auc_res_list)[i]
	title <- bquote(bold(.(paste0(LETTERS[i], ". "))) ~ "The proportion of causal EWAS genes = " ~ .(pec))
	plot_auc_res(auc_res_list[[i]], pec, 2, title_nam = title)	
})
names(auc_plot_list) <- unique_pec
glist <- lapply(auc_plot_list, ggplotGrob)

out_path <- "results/plots/method_test_gene_up_all"
system(paste("mkdir -p", out_path))

lapply(seq_along(auc_plot_list), function(x) {
	out_nam <- file.path(out_path, paste0("PEC_", names(auc_plot_list)[x], ".pdf"))
	p_out <- auc_plot_list[[x]]
	ggsave(out_nam, plot = p_out)
})

pdf("results/plots/methods_test_gene_up_auc_plot_all_databases.pdf")
marrangeGrob(glist, nrow=1, ncol=1, top = NULL)
dev.off()

## WORKOUT n_ewas_genes and n_gwas_genes
sim2_res <- sim2_res %>%
	mutate(n_gwas_genes = n_genes * gwas_power, 
		   n_ewas_genes = n_genes * ewas_power)


summary(sim2_res)

#### Testing best way to present summary of data

### just percent ewas causal at 1 and 0.05 at n=500 & n=10000
small_auc_res_list <- lapply(auc_res_list[c("0.05", "1")], function(x) {
	x %>%
		dplyr::filter(n_genes == 500 | n_genes == 10000)
})

i=1
small_auc_plot_list <- lapply(seq_along(small_auc_res_list), function(i) {
	pec <- names(small_auc_res_list)[i]
	title <- bquote(bold(.(paste0(LETTERS[i], ". "))) ~ "The proportion of causal EWAS genes = " ~ .(pec))
	plot_auc_res(small_auc_res_list[[i]], pec, 2, title_nam = title)	
})
glist <- lapply(small_auc_plot_list, ggplotGrob)
pdf("results/plots/methods_test_gene_up_auc_plot_all_databases_summary.pdf")
marrangeGrob(glist, nrow=2, ncol=1, top = NULL)
dev.off()

# -------------------------------------------------------
# Analyse the results for sim2 -- gene mapping compared to protein mapping
# -------------------------------------------------------

sim2_go_kegg_protein_res <- sim2_res_proteins %>%
	dplyr::filter(pathway_type %in% c("go", "kegg"))

sim2_res <- bind_rows(list(genes = sim2_res_genes, proteins = sim2_go_kegg_protein_res), 
					  .id = "genes_or_proteins") %>%
			mutate(pathway_type = paste0(pathway_type, "-", genes_or_proteins)) %>%
			dplyr::select(-genes_or_proteins)

unique_pdb <- unique(sim2_res$pathway_type)
variables <- c("or_g", "or_p", "rho_p")
colours <- brewer.pal(n = length(variables), name = "Set1")
unique_pec <- unique(sim2_res$percent_ewas_causal)
unique_pec <- unique_pec[unique_pec != 0]
# sim_res <- sim1_res

x=unique_pec[1]
auc_res_list <- lapply(unique_pec, function(x) {
	message("Proportion of EWAS hits that are causal: ", x)
	sim_res <- sim2_res %>%
		dplyr::filter(percent_ewas_causal == x | percent_ewas_causal == 0)
	suppressMessages(auc_res <- gen_auc_res(sim_res, 2))
})
names(auc_res_list) <- unique_pec
n_genes_labs <- paste("n_genes:", unique(sim2_res$n_genes))
names(n_genes_labs) <- as.character(unique(sim2_res$n_genes))

auc_plot_list <- lapply(seq_along(auc_res_list), function(i) {
	pec <- names(auc_res_list)[i]
	title <- bquote(bold(.(paste0(LETTERS[i], ". "))) ~ "The proportion of causal EWAS genes = " ~ .(pec))
	plot_auc_res(auc_res_list[[i]], pec, 2, title_nam = title)	
})
names(auc_plot_list) <- unique_pec

out_path <- "results/plots/method_test_gene_v_protein"
system(paste("mkdir", out_path))

lapply(seq_along(auc_plot_list), function(x) {
	out_nam <- file.path(out_path, paste0("PEC_", names(auc_plot_list)[x], ".pdf"))
	p_out <- auc_plot_list[[x]]
	ggsave(out_nam, plot = p_out)
})


glist <- lapply(auc_plot_list, ggplotGrob)
pdf("results/plots/methods_test_gene_up_auc_plot_go_and_kegg_all_genes.pdf")
marrangeGrob(glist, nrow=1, ncol=1, top = NULL)
dev.off()


# -------------------------------------------------------
# Checking the distribution of genes in each pathway for each dataset
# -------------------------------------------------------

go_terms <- read_tsv("data/derived/go_terms.txt") %>% 
	dplyr::filter(!is.na(pathway_id))

kegg_terms <- read_tsv("data/derived/kegg_terms.txt") 

stringdb_terms <- read_tsv("data/derived/stringdb_terms.txt")

reactome_terms <- read_tsv("data/derived/reactome_terms.txt")
ppi_terms <- read_tsv("data/derived/ppi_terms.txt")


get_n_genes <- function(database) 
{
	### get genes per pathway 
	pathway_data <- get(paste0(database, "_terms"))
	pathway_data %>%
		group_by(pathway_id) %>%
		summarise(n_genes = n())
}

databases <- c("kegg", "go", "stringdb", "reactome", "ppi")
n_gene_plots <- lapply(databases, function(db) {
	df <- get_n_genes(db)
	print(summary(df))
	ggplot(df, aes(x = n_genes)) + 
		geom_histogram(fill = "blue", colour = "black") +
		labs(title = db) +
		theme_bw()
})

glist <- lapply(n_gene_plots, ggplotGrob)
pdf("results/plots/genes_per_pathway_distributions.pdf")
marrangeGrob(glist, nrow=1, ncol=1, top = NULL)
dev.off()
