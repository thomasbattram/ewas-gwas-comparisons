# ------------------------------------------------------
# Extract string db data
# ------------------------------------------------------

# BiocManager::install("STRINGdb")
pkgs <- c("tidyverse", "STRINGdb")
lapply(pkgs, require, character.only = TRUE)

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

## NOTE: Did try using the stringdb package but there was an error
##		 and the package was a bit of a pain to use so ignoring that
##		 shite and just downloading the data

all_genes <- read_tsv("data/derived/ensembl_genes.txt", guess_max = 1e6)
prot_links <- data.table::fread("data/stringdb/9606.protein.links.full.v11.0.txt.gz")
prot_aliases <- data.table::fread("data/stringdb/9606.protein.aliases.v11.0.txt.gz")
colnames(prot_aliases) <- c("string_protein_id", "alias", "source")

# ------------------------------------------------------
# Extracting protein networks
# ------------------------------------------------------

str(prot_links)

prot_links %>%
	dplyr::filter(combined_score == 999) %>%
	head()

thresholds <- seq(399, 999, 100)
i <- thresholds[1]
n_protein_check <- lapply(thresholds, function(i) {
	df <- prot_links %>%
		dplyr::filter(combined_score > i)
	proteins <- c(df$protein1, df$protein2)
	uni_p <- unique(proteins)
	tab_p <- table(proteins)
	mean_p <- mean(tab_p)
	median_p <- median(tab_p)
	min_p <- min(tab_p)
	max_p <- max(tab_p)
	out_tab <- data.frame(threshold=i, 
						  proteins_left = length(uni_p), 
						  total_interactions = length(proteins), 
						  mean_interactions = mean_p, 
						  median_interactions = median_p, 
						  min_interactions = min_p, 
						  max_interactions = max_p)
	return(out_tab)
})
n_protein_check <- bind_rows(n_protein_check)
n_protein_check

# using threshold of 900+ 
threshold_to_use <- 899

prot_link_filtered <- prot_links %>%
	dplyr::filter(combined_score > threshold_to_use)

# ------------------------------------------------------
# Link proteins to ensembl genes
# ------------------------------------------------------

### Want a table like this
# pathway_id ensembl_gene_id description(optional)
##

unique(prot_aliases$source)[grep("ensembl_gene", unique(prot_aliases$source), ignore.case = T)]
unique(prot_aliases$source)[grep("hgnc", unique(prot_aliases$source), ignore.case = T)]
head(prot_aliases)
source2eb <- unlist(strsplit(unique(prot_aliases$source), " "))
head(source2eb)
unique(source2eb)

unique(prot_aliases$source)[grep("Ensembl_gene", unique(prot_aliases$source))]

stringdb_ensembl_genes <- prot_aliases %>%
	dplyr::filter(source == "Ensembl_ArrayExpress Ensembl_Source Ensembl_gene")

stringdb_ensembl_genes2 <- prot_aliases %>%
	dplyr::filter(source == "Ensembl_ArrayExpress Ensembl_HGNC_Ensembl_ID(supplied_by_Ensembl) Ensembl_Source Ensembl_gene")

sum(stringdb_ensembl_genes2$alias %in% stringdb_ensembl_genes$alias)
# 0
sum(stringdb_ensembl_genes2$alias %in% all_genes$ensembl_gene_id)
# 18373 (out of 18408)
sum(stringdb_ensembl_genes$alias %in% all_genes$ensembl_gene_id)
# 937 (out of 1158)
str(stringdb_ensembl_genes)

prot_aliases_filtered <- bind_rows(stringdb_ensembl_genes, stringdb_ensembl_genes2) %>%
	dplyr::select(-source) %>%
	distinct()

dim(prot_aliases_filtered)


tab1 <- table(prot_link_filtered$protein1)
tab2 <- table(prot_link_filtered$protein2)
head(tab1)
head(tab2)
# checking things are duplicated 
prot_link_filtered %>%
	dplyr::filter(protein1 == "9606.ENSP00000001146", protein2 == "9606.ENSP00000337915")
prot_link_filtered %>%
	dplyr::filter(protein1 == "9606.ENSP00000337915", protein2 == "9606.ENSP00000001146")
# They are!! --> good news as means can just remove a protein column

prot_out <- prot_link_filtered %>%
	dplyr::select(protein1) %>%
	left_join(prot_aliases_filtered, by = c("protein1" = "string_protein_id")) %>%
	dplyr::filter(alias %in% all_genes$ensembl_gene_id) %>%
	rename(pathway_id = protein1, ensembl_gene_id = alias)

# quick look at number of genes per pathway
prot_out %>%
	group_by(pathway_id) %>%
	summarise(n()) %>%
	summary
# maybe some outliers...

write.table(prot_out, file = "data/derived/stringdb_terms.txt", 
			col.names = T, row.names = F, quote = F, sep = "\t")

# FIN!













