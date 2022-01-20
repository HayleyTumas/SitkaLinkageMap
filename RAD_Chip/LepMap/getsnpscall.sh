cd #DataFiles


#Full QC SNPs (IND, GENO, and MAF)
cat RAD_SNP_RPCook4N_f.call|cut -f 1,2|awk '(NR>=7)' >RADSNPC4N_callfile.txt

#get SNPS from vcf to get SNP ID
grep -v "^##" RAD_SNP_2FAM_Cook4N.recode.vcf | cut -f1-3>RADSNPC4N_SNPIDs.txt
