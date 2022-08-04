# -----------------------------------
# dividing genome into chunks
# -----------------------------------

## This script splits the genome into chunks and assigns CpG sites 
## and SNPs to those chunks

pkgs <- c("tidyverse")
lapply(pkgs, require, character.only = T)

# -----------------------------------
# functions
# -----------------------------------

divide_chromosome <- function(data, start_bp, end_bp, window) 
{
	### splits chromosome into chunks based on the window size
	### and data entered
	bp_start <- start_bp - window/2
	if (sign(bp_start) == -1) bp_start <- 1
	divisions <- seq(bp_start, end_bp + window/2, window)
	out <- lapply(seq_along(divisions), function(section) {
		div_1 <- divisions[section]
		div_2 <- divisions[section + 1]
		if (section == length(divisions)) return(NULL)
		
		groups_only <- data.frame(group = section, positions = paste(div_1, div_2, sep = ":"), 
								  CHR = unique(data$CHR))
		out_dat <- data %>%
			dplyr::filter(bp < div_2 & bp >= div_1) %>%
			mutate(group = section) %>%
			mutate(positions = paste(div_1, div_2, sep = ":"))
		return(list(out_dat = out_dat, groups_only = groups_only))
	})
	out2 <- plyr::compact(out)
	out2 <- transpose(out2) %>% 
                 map(bind_rows)
	return(list(div_genome = out2$out_dat, groups_only = out2$groups_only))
}

# -----------------------------------
# Divide up the 450k data
# -----------------------------------

# fdata
load("~/tb_MB009_LungCancer/RDATA/fdata_new.RData")

# input data
# SNP/CPG name, CHR, bp

input_dat <- fdata.new %>%
	mutate(CpG = as.character(TargetID)) %>%
	mutate(bp = as.numeric(as.character(COORDINATE_37))) %>%
	dplyr::filter(!CHR %in% c("X", "Y")) %>%
	mutate(CHR = as.numeric(as.character(CHR))) %>%
	dplyr::select(CpG, CHR, bp) %>%
	dplyr::filter(!is.na(CHR))

chr <- unique(input_dat$CHR)[1]
chr <- 16
div_dat <- lapply(unique(input_dat$CHR), function(chr) {
	print(chr)
	dat <- input_dat %>%
		dplyr::filter(CHR == chr)

	out <- divide_chromosome(dat, start_bp = min(dat$bp), end_bp = max(dat$bp), window = 5e5)
	return(out)
})
div_dat <- transpose(div_dat) %>%
	map(bind_rows)
div_genome <- div_dat$div_genome
groups_dat <- div_dat$groups_only

# sort the group data
groups_dat <- groups_dat %>%
	arrange(CHR) %>%
	mutate(group_total = 1:nrow(.))

write.table(groups_dat, file = "data/derived/groups.txt",
			col.names = T, row.names = F, quote = F, sep = "\t")

# join the group data and CpG data so it now has "group_total" column
fin_div_genome <- div_genome %>%
	left_join(groups_dat)

write.table(fin_div_genome, file = "data/derived/epigenetic/divided_genome.txt",
			sep = "\t", row.names = F, col.names = T, quote = F)

# -----------------------------------
# Divide up the SNP data from 1k-genomes
# -----------------------------------

converter <- data.table::fread("~/AD/finemapping/rs_convert_2.txt")
snp_pos_dat <- converter %>%
	mutate(CHR = gsub(":.*", "", `1000G_MARKER`)) %>%
	dplyr::filter(!CHR %in% c("X", "Y")) %>%
	mutate(bp = gsub("^[0-9]*:", "", `1000G_MARKER`)) %>%
	mutate(bp = gsub(":.*$", "", bp)) %>%
	mutate(bp = as.numeric(bp)) 
rm(converter)
# write.table(snp_pos_dat, file = "data/snp_pos_converter_1KGen.txt", quote = F, row.names = F, col.names = T, sep = "\t")

nrow(dplyr::filter(snp_pos_dat, CHR == chr))
gen_dat <- map_dfr(unique(snp_pos_dat$CHR), function(chr) {
	print(chr)
	div_dat <- div_genome %>%
		dplyr::filter(CHR == chr)

	start_pos <- min(as.numeric(gsub(":.*", "", div_dat$positions)))
	end_pos <- max(as.numeric(gsub(".*:", "", div_dat$positions)))

	dat <- snp_pos_dat %>%
		dplyr::filter(CHR == chr) %>%
		dplyr::filter(bp <= end_pos & bp >= start_pos)

	out_dat <- divide_chromosome(dat, start_bp = min(div_dat$bp), end_bp = max(div_dat$bp), window = 5e5)
    return(out_dat$div_genome)
})

gen_dat <- gen_dat %>%
	mutate(CHR = as.numeric(CHR)) %>%
	left_join(groups_dat) %>%
	dplyr::filter(group_total %in% unique(fin_div_genome$group_total))

write.table(gen_dat, file = "data/derived/genetic/genetic_data_divided_genome.txt",
			sep = "\t", row.names = F, col.names = T, quote = F)
