# A comparison of EWAS and GWAS signal

## GitHub/lab notes

If pushing to private gitlab directory make sure to be in the "master" branch and push to "origin"
If pushing to public github directory make sure to be in the "public" branch and push to "github" like so 

`
git push github public:master
`

## Aims
1. Compare the regions identified by GWAS and EWAS
2. Assess if the genes/pathways identified by GWAS are also identified by EWAS more than expected by chance
3. Assess whether EWAS signal tends to correlate more with GWAS of non-corresponding traits

## Methods

To re-do the work and add more traits here are the steps you can skip: 1, 4, 5, 6, 7, 9, 10. Step 15 could potentially be skipped as that is extracting all data from the ieu open gwas project database so doesn't impact on data for selected traits. HOWEVER, it should be repeated as there would likely have been an update to the database.

If updating the data for a trait (e.g. if a new EWAS/GWAS of BMI comes out) will need to delete some of the data from various steps in the process!

### Section 1: Data extraction and setup
1. Divide genome up into 500kb blocks
	+ [divide_genome.R](R/divide_geneome.R)
2. Extract EWAS with N>4500
	+ [extract_ewas_data.R](R/extract_ewas_data.R) and [extract_ewas_data_manual.R](R/extract_ewas_data_manual.R)
3. Extract corresponding GWAS
	+ [extract_gwas_data.R](R/extract_gwas_data.R)
4. Extract genes and pathways
	+ [extract_genes_and_pathways.R](R/extract_genes_and_pathways.R)
5. Map SNP/CpG positions to genes and genes to pathways
	+ [mapping.R](R/mapping.R)
6. Extract data from epigraphdb 
	+ [extracting_epigraph_db_data.py](py/extracting_epigraph_db_data.py)
7. Clean up the epigraphdb data
	+ [clean_epigraphdb_data.R](R/clean_epigraphdb_data.R)

### Section 2: Direct overlap
8. For each trait, assess the overlap in GWAS and EWAS signals for each 500kb block
	+ [assess_overlap.R](R/assess_overlap.R)
9. Assess the gene overlap and check for commonalities across traits
	+ [gene_overlap.R](R/gene_overlap.R)

### Section 3: Gene/pathway overlap enrichment tests
10. Assess best method for conducting overlap tests and probable scenarios in which there is power to detect overlap
	+ [comparing_overlap_methods_sims.R](R/comparing_overlap_methods_sims.R) and [comparing_overlap_methods_sims_results.R](R/comparing_overlap_methods_sims_results.R)
11. Extract genes for empirical enrichment analyses
	+ [extract_genes_for_enrichment.R](R/extract_genes_for_enrichment.R)
12. Use best method to assess if genes/pathways tagged by GWAS are more likely to be identified by EWAS than expected by chance
	+ [enrichment_tests.R](R/enrichment_tests.R) and [enrichment_tests_results.R](R/enrichment_tests_results.R)
13. What are the different genesets that are enriched for? 
	+ [examining_enriched_pathways.R](R/examining_enriched_pathways.R)
14. Can we learn about the underlying shared architecture of traits from the overlap? 
	+ [architecture_sims.R](R/architecture_sims.R) and [architecture_sims_results.R](R/architecture_sims_results.R)

### Section 4: Gene/pathway overlap tests between all GWAS
15. Extract data from open gwas database
	+ [extract_opengwas_db_data.R](R/extract_opengwas_db_data.R)
16. Explore overlap between genes/pathways identified by GWAS signal from a large number of traits and the EWAS genes/pathways
	+ [overlapping_many_gwas_signals.R](R/overlapping_many_gwas_signals.R), [calulating_pathway_enrichment_correlations.R](R/calulating_pathway_enrichment_correlations.R) and [assessing_pathway_enrichment_correlations.R](R/assessing_pathway_enrichment_correlations.R)
