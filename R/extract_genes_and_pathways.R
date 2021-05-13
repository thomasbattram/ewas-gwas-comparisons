# ------------------------------------------------------------
# extracting genes and pathways
# ------------------------------------------------------------

## Aims of script:
## Extract all ensembl gene positions and link them to GO and KEGG
## terms

## pkgs
# library(tidyverse) # tidy code and data -- CAN'T USE WITH BIOMART FUNCTIONS!! -- CRASHES
library(FDb.InfiniumMethylation.hg19) # HM450 genomic annotations
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # HM450 genomic annotations
library(biomaRt) # ensembl gene IDs and positions
library(limma) # kegg pathways linked to genes

## directory
dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

# ------------------------------------------------------------
# extracting the genes, positions and pathways from biomart - test
# ------------------------------------------------------------
# mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# attri <- listAttributes(mart)
# chr_num <- "y"
# print(chr_num)
# ens_genes <- getBM(
# 	attributes=c("ensembl_gene_id","hgnc_symbol","uniprot_gn_id","chromosome_name","start_position","end_position"),
# 	filters=c("chromosome_name", "start", "end"),
# 	values=list(chromosome=chr_num, start="1",end=as.character(1e99)),
# 	mart=mart
# )
# ens_genes$chromosome_name <- as.character(ens_genes$chromosome_name)
# head(ens_genes)

# ens_genes[1,2]

# ids <- ens_genes %>%
# 	dplyr::filter(uniprot_gn_symbol != "")

# x <- ids$hgnc_symbol == ids$uniprot_gn_symbol
# ids[!x,]

# hgnc_genes <- getBM(
# 	attributes=c("hgnc_symbol","chromosome_name","start_position","end_position"),
# 	filters=c("chromosome_name", "start", "end"),
# 	values=list(chromosome=chr_num, start="1",end=as.character(1e99)),
# 	mart=mart
# )

# dim(hgnc_genes)
# dim(ens_genes)

# all(hgnc_genes$hgnc_symbol == ens_genes$hgnc_symbol) 
# Extraction works same way regardless

# ------------------------------------------------------------
# extracting genes from biomaRt
# ------------------------------------------------------------

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
attri <- listAttributes(mart)
chr <- c(1:22, "x", "y")
ori_time <- proc.time()
all_genes <- lapply(chr, function(chr_num) {
	print(chr_num)
	genes <- getBM(
	    attributes=c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position"),
	    filters=c("chromosome_name", "start", "end"),
	    values=list(chromosome=chr_num, start="1",end=as.character(1e99)),
    	mart=mart
    )
    genes$chromosome_name <- as.character(genes$chromosome_name)
	return(genes)   
})
library(tidyverse) # couldn't load this before as it crashed things! 
all_genes <- bind_rows(all_genes) %>%
	distinct()
fin_time <- proc.time() - ori_time 
fin_time # elapsed = 7.949

write.table(all_genes, file = "data/derived/ensembl_genes.txt", 
			row.names = F, col.names = T, quote = F, sep = "\t")

# ------------------------------------------------------------
# extracting GO terms from biomaRt
# ------------------------------------------------------------

## go attributes
# go_id = the ID 
# name_1006 = short description of id
# definition_1006 = long description of id
# go_linkage_type = codes that describe the type of evidence used to assign
#					genes to certain go terms. For example IDA = Inferred from Direct Assay.
#					Guide to codes can be found here: http://geneontology.org/docs/guide-go-evidence-codes/
# namespace_1003 = the 3 go term domains:
#				   cellular component, biological process, molecular function

go_pathways <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006"),
                 mart = mart, verbose = T)
go_pathways <- go_pathways %>%
	dplyr::rename(description = name_1006) %>%
	dplyr::rename(pathway_id = go_id)

write.table(go_pathways, file = "data/derived/go_terms.txt", 
			row.names = F, col.names = T, quote = F, sep = "\t")

# ------------------------------------------------------------
# extracting KEGG terms from limma
# ------------------------------------------------------------

# extract kegg pathways from limma
kegg_pathways <- getGeneKEGGLinks(species.KEGG = "hsa", convert = TRUE)
# change gene IDs to ensembl_gene_ids 
gene_ids <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"), mart = mart)
head(gene_ids)
sum(unique(kegg_pathways$GeneID) %in% gene_ids$entrezgene_id) # missing ~170 (fine with this...)

# get kegg pathway descriptions
kegg_pathway_names <- getKEGGPathwayNames(species.KEGG = "hsa", remove.qualifier = TRUE)
head(kegg_pathway_names)

require(tidyverse)
# put it all together and write it out
kegg_pathways2 <- kegg_pathways %>%
	mutate(entrezgene_id = as.numeric(GeneID)) %>%
	left_join(gene_ids) %>%
	left_join(kegg_pathway_names) %>%
	dplyr::select(-entrezgene_id, -GeneID) %>%
	dplyr::rename(pathway_id = PathwayID) %>%
	dplyr::rename(description = Description)

write.table(kegg_pathways2, file = "data/derived/kegg_terms.txt", 
			row.names = F, col.names = T, quote = F, sep = "\t")
