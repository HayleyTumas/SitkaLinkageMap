cd #DataFiles


#Full QC SNPs (IND, GENO, and MAF)
cat RADQC_2FAMSSf.call|cut -f 1,2|awk '(NR>=7)' >2FAMSS_callfile.txt

#get SNPS from vcf to get SNP ID
grep -v "^##" RAD_comb2_QC_BestSNP.recode.vcf | cut -f1-3>2FAMSS_SNPIDs.txt
