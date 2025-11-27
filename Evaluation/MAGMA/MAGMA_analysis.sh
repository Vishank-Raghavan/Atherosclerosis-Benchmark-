#!/usr/bin/env bash

#Map SNPs in reference data to genes.
#Assumes you have added magma to your path
magma --annotate --snp-loc g1000_eur.bim --gene-loc NCBI38.gene.loc --out annotation

#Run Gene Analysis via MAGMA
magma --bfile g1000_eur \
        --pval GWAS_Summary_Stats.tsv use=hm_rsid,p_value N=184305 \
        --gene-annot annotation.genes.annot \
        --out coronary_artery_disease


#Filter genes for bonferroni corrected significance
awk 'NR==1 || $9 < 0.0000025' coronary_artery_disease.genes.out > significant_genes.txt


#If list is too small can also check using FDR
./fdr_adjustment.py
