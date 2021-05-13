

calc_increments <- function(length_task, percent_increments=10) 
{
	### Calculate increments to output message
	increments <- round(length_task / percent_increments)
	message_at <- round(seq(increments, length_task, length.out = percent_increments))
	names(message_at) <- seq(percent_increments, 100, percent_increments)
	return(message_at)
}

output_percent_complete <- function(n_task, increments) 
{
	### output percentage complete
	if (n_task %in% increments) {
		percent <- increments[increments == n_task]
		message(names(percent), "% complete.")
	}
}


generate_genesets <- function(all_pathways, pathways_and_genes, n_causal_pathways, n_consequent_pathways) {
	ca_pathways <- sample(x = all_pathways, size = n_causal_pathways, replace = FALSE)
	ca_genes <- unique(subset(pathways_and_genes, pathway_id %in% ca_pathways)$ensembl_gene_id)
	con_pathways <- sample(x = all_pathways, size = n_consequent_pathways, replace = FALSE)
	con_genes <- unique(subset(pathways_and_genes, pathway_id %in% con_pathways)$ensembl_gene_id)
	# message("number of causal genes: ", length(ca_genes))
	# message("number of consequent genes: ", length(con_genes))
	out <- list(causal = ca_genes, consequential = con_genes)
	return(out)
}

# overlap test function
overlap_test <- function(group1, group2, all_variables) 
{
    # asks if there is more of an overlap between
    # group one and group two than expected by chance if 
    # those groups were randomly sampled from all the variables

    q <- sum(group1 %in% group2) # overlap between ewas and gwas pathways
    m <- length(group2) - q # number of pathways identified by GWAS but not EWAS
    k <- length(group1) - q # number of pathways identified by EWAS but not by GWAS
    n <- length(all_variables) - q - m - k # all pathways minus pathways identified by EWAS and GWAS
    tab <- matrix(c(q, m, k, n), 2, 2)
    out <- fisher.test(tab, alternative = "greater")
    return(out)
}

extract_pathways <- function(gwas_genes, ewas_genes, databases)
{
	# function to extact pathways from each of the databases
	# supplied. Database dfs should be name DATABASE_terms

	out_paths <- lapply(seq_along(databases), function(x) {
		pathway_dat <- get(paste0(databases[x], "_terms"))
		gwasp <- pathway_dat %>%
			dplyr::filter(ensembl_gene_id %in% gwas_genes) %>%
			pull(pathway_id) %>%
			unique
		ewasp <- pathway_dat %>%
			dplyr::filter(ensembl_gene_id %in% ewas_genes) %>%
			pull(pathway_id) %>%
			unique
		return(list(gwas = gwasp, ewas = ewasp))
	})
	out_paths <- flatten(out_paths)
	return(out_paths)
}

extract_genes_in_pathways <- function(pathway_df) 
{
	# takes all genes from each pathway and outputs this as a list    

    unique_pathways <- unique(pathway_df$pathway_id)
    len_paths <- length(unique_pathways)
    all_pathways <- lapply(seq_along(unique_pathways), function(x) {
        message("Pathway ", x, " of ", len_paths)
    	pathway_of_interest <- unique_pathways[x]
    	pathway_genes <- pathway_df %>%
    		dplyr::filter(pathway_id %in% pathway_of_interest) %>%
    		pull(ensembl_gene_id)
    	return(pathway_genes)
    })
    names(all_pathways) <- unique_pathways
    return(all_pathways)
}

perform_enrichment <- function(identified_genes, database, background_genes) 
{
    # uses fishers exact test to perform geneset enrichment analyses 
    # for each pathway and output pathway_id and "enrichment scores"

	pathway_gene_data <- get(paste0(database, "_genes"))
	len_paths <- length(pathway_gene_data)
    enriched_paths <- map_dfr(seq_along(pathway_gene_data), function(x) {
        # message("Pathway ", x, " of ", len_paths)
        pathway_of_interest <- names(pathway_gene_data)[x]      
        pathway_genes <- pathway_gene_data[[pathway_of_interest]]

        res <- overlap_test(group1 = identified_genes, 
                            group2 = pathway_genes, 
                            all_variables = background_genes)
        out <- tibble::tibble(pathway_id = pathway_of_interest, enrich_or = res$estimate, 
        				  enrich_p = res$p) %>%
        		mutate(adj_enrich_p = p.adjust(enrich_p, method = "BH"))
        return(out)
    }) 
}

sim1 <- function(prop_causal_pathways, gwas_power, ewas_power, percent_ewas_causal, prop_consequent_pathways, pathway_type, sim)
{
	# extract all pathways and number of pathways
	set.seed(sim)
	all_terms <- get(paste0("u_", pathway_type)) 
	pathway_dat <- get(paste0(pathway_type, "_terms"))
	n_causal_pathways <- round(prop_causal_pathways * length(all_terms))
	n_consequent_pathways <- round(prop_consequent_pathways * length(all_terms))
	# extract genes from pathways
	genes <- generate_genesets(all_terms, pathway_dat, n_causal_pathways, n_consequent_pathways)
	# sample from the genes and assess gene overlap
	gwasg <- sample(genes$causal, gwas_power * length(genes$causal), replace=FALSE)
	ewasg <- sample(genes$causal, ewas_power * length(genes$causal) * percent_ewas_causal, replace=FALSE)
	ewasg <- unique(c(ewasg, sample(genes$consequential, ewas_power * length(genes$consequential) * (1 - percent_ewas_causal), replace=FALSE)))
	outg <- overlap_test(ewasg, gwasg, all_unique_genes) # gene overlap
	# map the genes back to pathways and assess pathway overlap
	gwasp <- unique(pathway_dat[pathway_dat$ensembl_gene_id %in% gwasg, "pathway_id", drop = T])
	ewasp <- unique(pathway_dat[pathway_dat$ensembl_gene_id %in% ewasg, "pathway_id", drop = T])
	outp <- overlap_test(ewasp, gwasp, all_terms) # pathway overlap
	# do pathway enrichment and correlation
	gwas_enrich <- perform_enrichment(gwasg, database = pathway_type, all_unique_genes)
	ewas_enrich <- perform_enrichment(ewasg, database = pathway_type, all_unique_genes)
	enrich_cor <- cor(log(gwas_enrich$enrich_or), log(ewas_enrich$enrich_or),
					  use = "pairwise.complete.obs",
					  method = "spearman")
	# write out results
	out <- list(gene = outg, 
				pathway = outp, 
				enrich_cor = enrich_cor, 
				n_gwas_genes = length(gwasg), 
				n_ewas_genes = length(ewasg), 
				gene_overlap = sum(ewasg %in% gwasg),
				pathway_overlap = sum(ewasp %in% gwasp))
	return(out)
}

apply_sim1 <- function(params) {
	out <- lapply(split(params, 1:nrow(params)), function(x) {
		# output_percent_complete(rownames(x), message_at)
		out <- sim1(x$prop_causal_pathways, x$gwas_power, x$ewas_power, 
					x$percent_ewas_causal, x$prop_consequent_pathways, x$pathway_type,
					x$sim)
		x$or_g <- out$gene$estimate
		x$or_p <- out$pathway$estimate
		x$rho_p <- out$enrich_cor
		x$n_gwas_genes <- out$n_gwas_genes
		x$n_ewas_genes <- out$n_ewas_genes
		x$gene_overlap <- out$gene_overlap
		x$pathway_overlap <- out$pathway_overlap
		return(x)
	})
	return(out)
}

sim2 <- function(n_genes, ewas_power, gwas_power, percent_ewas_causal, database, sim) 
{
	# simulation function

	set.seed(sim)
	pathway_dat <- get(paste0(database, "_terms"))
	# sample genes
	causal_trait_genes <- sample(all_unique_genes, n_genes)
	gwas_genes <- sample(causal_trait_genes, n_genes * gwas_power)
	n_ewas_genes <- n_genes * ewas_power
	n_ewas_causal_genes <- n_ewas_genes * percent_ewas_causal
	ewas_genes <- unique(c(sample(causal_trait_genes, n_ewas_causal_genes), 
					sample(all_unique_genes, n_ewas_genes - n_ewas_causal_genes)))
	# perform gene overlap test
	message("Testing if genes overlap")
	overlap_g <- overlap_test(ewas_genes, gwas_genes, all_unique_genes)
	# extract pathways and do overlap test on them all
	all_pathways <- extract_pathways(gwas_genes, ewas_genes, database)
	gwasp <- all_pathways$gwas
	ewasp <- all_pathways$ewas
	message("Testing if pathways overlap")
	overlap_p <- overlap_test(ewasp, gwasp, unique(pathway_dat$pathway_id))
	# perform enrichment analyses
	message("Performing enrichment analyses")
	gwas_enrich <- perform_enrichment(gwas_genes, database, all_unique_genes)
	ewas_enrich <- perform_enrichment(ewas_genes, database, all_unique_genes)
	enrich_cor <- cor(log(gwas_enrich$enrich_or), log(ewas_enrich$enrich_or),
					  use = "pairwise.complete.obs",
					  method = "spearman")

	out <- list(gene_overlap = sum(ewas_genes %in% gwas_genes), 
				gene_res = overlap_g$estimate, 
				pathway_overlap = sum(ewasp %in% gwasp), 
				or_p = overlap_p$estimate, 
				rho_p = enrich_cor)
	return(out)
}

apply_sim2 <- function(params) {
	out <- lapply(split(params, 1:nrow(params)), function(x) {
		output_percent_complete(rownames(x), message_at)
		sim_res <- sim2(n_genes = x$n_genes, 
						ewas_power = x$ewas_power, 
						gwas_power = x$gwas_power, 
						percent_ewas_causal = x$percent_ewas_causal, 
						database = x$pathway_type, 
						sim = x$sim)
		x$gene_overlap <- sim_res$gene_overlap
		x$pathway_overlap <- sim_res$pathway_overlap
		x$or_g <- sim_res$gene_res
		x$or_p <- sim_res$or_p
		x$rho_p <- sim_res$rho_p
		return(x)
	})
	return(out)
}



sim3 <- function(prop_causal_genes, n_ewas_genes, n_gwas_genes, percent_ewas_causal, prop_consequent_genes)
{
	# 1. sample causal and consequential genes from total number of genes
	# 2. sample gwas and ewas genes and do gene overlap tests
	# 3. link these to pathways and do pathway overlap tests
	# 4. write out results

	# 1.
	ca_genes <- sample_n(all_genes, prop_causal_genes * nrow(all_genes))
	con_genes <- sample_n(all_genes, prop_consequent_genes * nrow(all_genes))

	# 2. 
	gwasg <- sample(ca_genes$ensembl_gene_id, n_gwas_genes)
	ewas_cag <- sample(ca_genes$ensembl_gene_id, n_ewas_genes * percent_ewas_causal)
	ewasg <- unique(c(ewas_cag, sample(con_genes$ensembl_gene_id, n_ewas_genes - length(ewas_cag))))

	outg <- overlap_test(ewasg, gwasg, all_genes$ensembl_gene_id)

	# 3.
	database <- c("go", "kegg") 
	outp <- lapply(1:2, function(x) {
		pathway_dat <- get(paste0(database[x], "_terms"))
		gwasp <- pathway_dat %>%
			dplyr::filter(ensembl_gene_id %in% gwasg) %>%
			pull(pathway_id) %>%
			unique
		ewasp <- pathway_dat %>%
			dplyr::filter(ensembl_gene_id %in% ewasg) %>%
			pull(pathway_id) %>%
			unique
		out_res <- overlap_test(gwasp, ewasp, unique(pathway_dat$pathway_id))
		overlap <- sum(ewasp %in% gwasp)
		return(list(res = out_res, overlap = overlap))
	})
	names(outp) <- database
	
	# 4.
	out <- list(gene = outg, 
				gene_overlap = sum(ewasg %in% gwasg),
				go = outp$go$res,  
				go_overlap = outp$go$overlap, 
				kegg = outp$kegg$res, 
				kegg_overlap = outp$kegg$overlap)
	return(out)
}

apply_sim3 <- function(params) {
	out <- lapply(split(params, 1:nrow(params)), function(x) {
		out <- sim3(x$prop_causal_genes, x$n_ewas_genes, x$n_gwas_genes, 
					x$percent_ewas_causal, x$prop_consequent_genes)
		x$gene_overlap <- out$gene_overlap
		x$go_overlap <- out$go_overlap
		x$kegg_overlap <- out$kegg_overlap
		x$or_g <- out$gene$estimate
		x$or_go <- out$go$estimate
		x$or_kegg <- out$kegg$estimate
		x$p_g <- out$gene$p.value
		x$p_go <- out$go$p.value 
		x$p_kegg <- out$kegg$p.value
		return(x)
	})
	return(out)
}


# n_genes <- total_genes 
# ewas_genes <- unique(gene_dat$ewas_genes)
# gwas_genes <- unique(gene_dat$gwas_genes)
# ca_con_overlap <- x$ca_con_overlap

sim4 <- function(n_genes, ewas_genes, gwas_genes, ca_con_overlap)
{
	# 1. get causal genes from gwas genes and sample rest from all genes
	# 2. get consequential genes from ewas genes and sample rest from 
	#	 causal genes * ca_con_overlap + all_genes - (causal genes * ca_con_overlap)
	# 3. Perform gene overlap tests
	# 4. Perform pathway overlap tests
	# 5. Write out results

# N genes (NG)
# gwas genes (GG)
# ewas genes (EG)
# overlapping genes (OG)
# gwas only genes (GOG)
# ewas only genes (EOG)
# casual ewas genes (CaEG)
# consequential ewas genes (CaEG)
# leftover genes (LG)

# need to work out:
# causal genes (CaG)
# consequential genes (ConG) 

# vary overlap 

# GOG = GG - OG
# EOG = EG - OG

# CaEG = OG + sample(EOG, l(EOG) * overlap)
# CaG = CaEG + GOG + sample(LG, NG - l(CaEG) - l(GOG))

# ConEG = EG - CaEG
# ConCaG = sample(CaG - CaEG, l(.) * overlap) + CaEG
# ConG = ConCaG + ConEG + sample(LG - CaG - ConEG, NG - l(ConCaG) - l(ConEG))



# ca genes = gwas genes + sample(non_ge_genes, n)

	# 0. 

	overlapping_g <- ewas_genes[ewas_genes %in% gwas_genes]
	ewas_only_g <- ewas_genes[!ewas_genes %in% overlapping_g]
	gwas_only_g <- gwas_genes[!gwas_genes %in% overlapping_g]
	ewas_ca_genes <- c(sample(ewas_only_g, length(ewas_only_g) * ca_con_overlap),
					   overlapping_g)
	ewas_con_genes <- ewas_genes[!ewas_genes %in% ewas_ca_genes]
	
	non_ge_genes <- all_genes %>%
		dplyr::filter(!ensembl_gene_id %in% ewas_genes &
					  !ensembl_gene_id %in% gwas_genes)


	# 1.
	n_gwas_genes <- length(gwas_genes)
	gwas_ewas_gene_dat <- all_genes %>%
		dplyr::filter(ensembl_gene_id %in% c(ewas_ca_genes, gwas_only_g))

	ca_genes <- non_ge_genes %>%
			sample_n(n_genes - nrow(gwas_ewas_gene_dat))
	ca_genes <- bind_rows(list(ca_genes, gwas_ewas_gene_dat))

	# 2.
	n_ewas_genes <- length(ewas_genes)

	max_overlap <- max(c(n_genes * ca_con_overlap, length(ewas_ca_genes)))

	con_ca_g <- ca_genes %>%
		dplyr::filter(!ensembl_gene_id %in% ewas_ca_genes) %>%
		sample_n(max_overlap - length(ewas_ca_genes)) %>%
		bind_rows(dplyr::filter(ca_genes, ensembl_gene_id %in% ewas_ca_genes))

	con_eg <- all_genes %>%
		dplyr::filter(ensembl_gene_id %in% ewas_con_genes)

	leftover_con_genes <- non_ge_genes %>%
		dplyr::filter(!ensembl_gene_id %in% c(ca_genes$ensembl_gene_id, 
											  con_eg$ensembl_gene_id)) %>%
		sample_n(n_genes - nrow(con_ca_g) - nrow(con_eg))

	con_genes <- bind_rows(list(con_ca_g, con_eg, leftover_con_genes))

	ewas_g <- sample(con_genes$ensembl_gene_id, n_ewas_genes)
	gwas_g <- sample(ca_genes$ensembl_gene_id, n_gwas_genes)

	# 3. 
	outg <- overlap_test(ewas_g, gwas_g, all_genes$ensembl_gene_id)

	# 4a.
	database <- "go"
	all_pathways <- extract_pathways(gwas_g, ewas_g, database)
	gwasp <- all_pathways$gwas
	ewasp <- all_pathways$ewas
	overlap_p <- overlap_test(ewasp, gwasp, unique(go_terms$pathway_id))
	# 4b.
	gwas_enrich <- perform_enrichment(gwas_g, database, all_genes$ensembl_gene_id)
	ewas_enrich <- perform_enrichment(ewas_g, database, all_genes$ensembl_gene_id)
	enrich_cor <- cor(log(gwas_enrich$enrich_or), log(ewas_enrich$enrich_or),
					  use = "pairwise.complete.obs",
					  method = "spearman")

	# 5.
	out <- list(gene = outg, 
				gene_overlap = length(overlapping_g),
				pathway_overlap = sum(ewasp %in% gwasp),
				or_g = outg$estimate,
				or_p = overlap_p$estimate,
				rho_p = enrich_cor)
	return(out)
}


# x <- params[1,]
apply_sim4 <- function(params) {
	out <- lapply(split(params, 1:nrow(params)), function(x) {
		print(x)
		x$trait <- as.character(x$trait)
		gene_dat <- trait_gene_dat %>%
			dplyr::filter(trait == x$trait)
		n_gene_dat <- trait_n_gene_dat %>%
			dplyr::filter(trait == x$trait)
		total_genes <- n_gene_dat$n_ewas_genes * x$n_gene_proportions +
					   n_gene_dat$n_gwas_genes * x$n_gene_proportions
		if (total_genes > length(all_genes$ensembl_gene_id)/2) return(NULL)
		out <- sim4(total_genes, 
					unique(gene_dat$ewas_genes), unique(gene_dat$gwas_genes),
					x$ca_con_overlap)
		x$gene_overlap <- out$gene_overlap
		x$or_g <- out$or_g
		x$or_p <- out$or_p
		x$rho_p <- out$rho_p
		return(x)
	})
	return(out)
}

