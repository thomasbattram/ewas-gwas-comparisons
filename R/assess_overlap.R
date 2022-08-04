# -------------------------------------------------
# Assessment of overlap between EWAS and GWAS
# -------------------------------------------------

## This script takes the mapped CpG sites and SNPs and looks 
## at the overlap in identified DMPs and SNPs from EWAS and GWAS
## with regards to genomic position

pkgs <- c("tidyverse", "pROC", "RColorBrewer", "ggrepel", "usefunc")
lapply(pkgs, require, character.only = T)

home_dir <- "~/projects/epi_gen_comp" # CHANGE ME WHEN NEEDED
source(file.path(home_dir, "R/mapping_functions.R"))

traits <- readLines("data/traits.txt")
check_all_dirs(traits)

# read in the data
# bmi_dat <- read_derived_data(trait = "body_mass_index")
# edu_dat <- read_derived_data(trait = "educational_attainment")
# alc_dat <- read_derived_data(trait = "alcohol_consumption_per_day")

cpg_divisions <- read_tsv("data/derived/epigenetic/divided_genome.txt", guess_max = 1e6)


# -------------------------------------------------
# does GWAS p predict EWAS hit? 
# -------------------------------------------------
generate_hit_dat <- function(gen_dat, meth_dat) {
	gen_dat_stats <- gen_dat %>%
		group_by(group_total) %>%
		summarise(min_p = min(p), med_p = median(p), 
				  max_beta = max(b), med_beta = median(b))

	dnam_dat_hits <- data.frame(group_total = 1:max(cpg_divisions$group_total)) %>%
		mutate(hit = ifelse(group_total %in% meth_dat$group_total, 1, 0))

	comb_dat <- gen_dat_stats %>%
		left_join(dnam_dat_hits)
	return(comb_dat)
}

get_auc <- function(gen_dat, meth_dat, out_roc = FALSE) {
	roc_all_dat <- generate_hit_dat(gen_dat, meth_dat)
	if (!1 %in% roc_all_dat$hit) return(NULL)
	auc_out <- lapply(seq_along(roc_all_dat)[-c(1,ncol(roc_all_dat))], function(x) {
		predictor <- colnames(roc_all_dat)[x]
		fom <- as.formula(paste("hit ~", predictor))
		roc_res <- roc(fom, data = roc_all_dat)
		if (out_roc) return(roc_res)

		auc <- t(as.data.frame(ci.auc(roc_res)))
		auc <- as.data.frame(auc) %>%
			mutate(trait = unique(meth_dat$Trait)) %>%
			mutate(predictor = predictor)
		colnames(auc)[1:3] <- c("auc_ci_low", "auc", "auc_ci_upper")
		return(auc)
	})
	if (!out_roc) auc_out <- do.call(rbind, auc_out)
	return(auc_out)

}
x <- 2
auc_dat <- lapply(traits, function(trait) {
	print(trait)
	dat <- read_derived_data(trait = trait)
	auc <- get_auc(gen_dat = dat$gen, meth_dat = dat$meth)
	return(auc)
})

### NEED TO NOTE INSULIN IS MISSING BECAUSE NO OVERLAP IN REGIONS!

fin_auc_dat <- do.call(rbind, auc_dat)
write.table(fin_auc_dat, file.path(home_dir, "report/report_data/auc_data.txt"), row.names = F, col.names = T, quote = F, sep = "\t")


# b_p <- c("min_p", "max_beta")
# auc_plot <- ggplot(dplyr::filter(fin_auc_dat, predictor %in% b_p), aes(x = trait, y = auc, colour = predictor)) +
#   geom_point(position = position_dodge(width = 0.9)) +
#   geom_linerange(aes(ymin = auc_ci_low, ymax = auc_ci_upper), position = position_dodge(width = 0.9)) +
#   labs(x = bquote("Trait"), 
#        y = "Area under the curve") +
#   scale_y_continuous(limits = c(0.4, 1))

# ggsave("results/temp/auc_scatter_minp.pdf", auc_plot)

# -------------------------------------------------
# What is the physical overlap?
# -------------------------------------------------
# new barplot res! 
bar_res <- map_dfr(traits, function(trait) {
	print(trait)
	# load trait
	dat <- read_derived_data(trait = trait)

	hits_df <- data.frame(group = 1:max(cpg_divisions$group_total), epi_g = NA, gen_g = NA)
	pval_threshold <- ifelse(max(dat$meth$P) < 1e-5, max(dat$meth$P), 1e-5)
	epi_group <- dat$meth %>%
		dplyr::filter(P < pval_threshold)
	epi_group <- unique(epi_group$group_total)
	gen_group <- dat$gen %>%
		dplyr::filter(p < pval_threshold)
	gen_group <- unique(gen_group$group_total)

	hits_df <- map_dfr(1:nrow(hits_df), function(i) {
		df <- hits_df[i, ] %>%
			mutate(epi_g = ifelse(i %in% epi_group, TRUE, FALSE), 
				   gen_g = ifelse(i %in% gen_group, TRUE, FALSE))
		return(df)
	})

	hits_df <- hits_df %>%
		mutate(trait = trait, sum_g = epi_g + gen_g) %>% 
		mutate(hits = case_when(epi_g == TRUE & gen_g == FALSE ~ "EWAS",
								epi_g == FALSE & gen_g == TRUE ~ "GWAS", 
								epi_g == TRUE & gen_g == TRUE ~ "both", 
								epi_g == FALSE & gen_g == FALSE ~ "neither")) %>%
		dplyr::filter(group %in% cpg_divisions$group_total)
	return(hits_df)
})

write.table(bar_res, file = "data/derived/data_for_barchart.txt", 
			col.names = T, row.names = F, quote = F, sep = "\t")

bar_res <- read_tsv("data/derived/data_for_barchart.txt")

# plotting it! 
bar_tab <- bar_res %>%
	group_by(trait, hits) %>%
	summarise(count = n()) %>%
	mutate(label_ypos = cumsum(count) - 0.5 * count) %>%
	arrange(trait, hits)

lapply(traits, function(tr) {
	all_hits <- bar_tab %>%
		dplyr::filter(trait == tr) %>%
		dplyr::filter(hits != "neither") %>%
		pull(count) %>%
		sum
	overlapping_hits <- bar_tab %>%
		dplyr::filter(trait == tr) %>%
		dplyr::filter(hits == "both") %>%
		pull(count)
	prop_overlapping <- overlapping_hits / all_hits
	return(prop_overlapping)
})

trait_labels <- c("alcohol_consumption_per_day" = "AC",
				  "body_mass_index" = "BMI", 
				  "educational_attainment" = "EA", 
				  "c-reactive_protein" = "CRP", 
				  "former_versus_never_smoking" = "FsNs", 
				  "current_versus_never_smoking" = "CsNs", 
				  "egfr" = "eGFR",
				  "urate" = "urate", 
				  "leisure_time_physical_activity" = "PA",  
				  "television_viewing_time" = "TV",				  
				  "fasting_glucose" = "Gluc", 
				  "fasting_insulin" = "Ins", 
				  "systolic_blood_pressure" = "SBP", 
				  "diastolic_blood_pressure" = "DBP", 
				  "birthweight" = "BW", 
				  "cognitive_abilities:_digit_test" = "Cog", 
				  "fev1" = "FEV1")

cols <- get_cb_palette()
cols <- cols[c(2:5)]
names(cols) <- c("both", "neither", "GWAS", "EWAS")
bar_p <- ggplot(bar_tab, aes(x = trait, y = count, fill = forcats::fct_rev(hits))) +
	geom_bar(stat = "identity") +
	# geom_text(aes(y = label_ypos, label = count), 
	# 		  colour = "black", size = 3) +
	scale_x_discrete("Trait", labels = trait_labels) +
	scale_fill_manual(values = cols) + 
	labs(y = "Regions") +
	theme_bw() +
	theme(legend.title = element_blank(), text = element_text(size = 13))

ggsave("results/plots/all_traits_overlap_bar.pdf", plot = bar_p, width = 10)
# should move manually to report section...
ggsave(file.path(home_dir, "report/report_data/figure_dat/all_traits_overlap_bar.pdf"), plot = bar_p, width=10)
