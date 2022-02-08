#!/bin/bash

YOURDIR=#DataFiles
cd $YOURDIR

##STEP 1: Merge separate family files into combined file
conda activate BCF
#this environment contains BCFtools

bcftools view fam1_qcedSNPs.recode.vcf -Oz -o FAM1_QC.vcf.gz
bcftools index FAM1_QC.vcf.gz

bcftools view fam2_qcedSNPs.recode.vcf -Oz -o FAM2_QC.vcf.gz
bcftools index FAM2_QC.vcf.gz

bcftools merge FAM1_QC.vcf.gz FAM2_QC.vcf.gz -Ov -o RAD_comb2_QC.vcf

conda deactivate

##STEP 2: Select Single 'Best SNP' by determining read depth and creating file with SNPs of highest depth
conda activate VCFTOOLS
#this environment contains vcftools

vcftools --vcf RAD_comb2_QC.vcf --site-depth --out Comb2_SNPS
#read in R to determined highest read depth per locus

#create new file with a single SNP per locus based on read depth
vcftools --vcf RAD_comb2_QC.vcf --snps BEST_COMB2_SNP.txt --out RAD_comb2_QC_BestSNP --recode
