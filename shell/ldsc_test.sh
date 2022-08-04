#!/bin/bash

## testing ldsc

wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1kg_eur.tar.bz2
tar -jxvf 1kg_eur.tar.bz2
ldsc.py \
--bfile 22 \
--l2 \
--ld-wind-cm 1 \
--out 22


### munging sum stats for smoking init

munge_sumstats.py \
--sumstats pgc.cross.SCZ17.2013-05.txt \
--N 17115 \
--out scz \
--merge-alleles w_hm3.snplist


munge_sumstats.py \
--sumstats ../gwas_sum_stats/liu_smoking_initiation_2019_sumstats.txt \
--snp RSID \
--N-col N \
--a1 ALT \
--a2 REF \
--p PVALUE \
--frq AF \
--ignore CHROM,POS,STAT,SE,EFFECTIVE_N,Number_of_Studies,ANNO,ANNOFULL \
--out smok_init \
--merge-alleles w_hm3.snplist

