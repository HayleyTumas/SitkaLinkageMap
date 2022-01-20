cd #DataFiles

conda activate VCFTOOLS

vcftools --vcf RAD_SNP_2FAMcomb.vcf --snps RADSNP_SCLM2_RemovedProb_Cook4n.txt  --out RAD_SNP_2FAM_Cook4N --recode
