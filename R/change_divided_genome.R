#######

dir <- "~/main_project/epi_gen_comp/"
setwd(dir)

library(tidyverse)

div_dat <- data.table::fread("/newhome/tb13101/main_project/epi_gen_comp/data/derived/genetic/genetic_data_divided_genome.txt")
print(pryr::object_size(div_dat))
div_dat <- div_dat %>%
	dplyr::select(-`1000G_MARKER`)
print(pryr::object_size(div_dat))
write.table(div_dat, file="data/derived/genetic/genetic_data_divided_genome.txt",
			row.names = F, col.names = T, quote = F, sep = "\t")
# save a subset so it's easier to work with! 
div_dat_sub <- div_dat[sample(1:nrow(div_dat), 1e6)]
print(table(div_dat_sub$CHR))

write.table(div_dat_sub, file="data/derived/genetic/sample_of_genetic_data_divided_genome.txt", 
			row.names = F, col.names = T, quote = F, sep = "\t")
