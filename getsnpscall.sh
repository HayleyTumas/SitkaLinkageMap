YOURDIR=/mnt/c/Users/dops0764/Documents/LepMap/Chip_QC/LM1
cd $YOURDIR

#Full QC SNPs (IND, GENO, and MAF)
cat CombQC_NRf.call|cut -f 1,2|awk '(NR>=7)' >CombQC_NR_snps_callfile.txt

#get SNPS from vcf to get SNP ID
grep -v "^##" Comb_QC_Noreps.recode.vcf | cut -f1-3>CombQC_NoReps_SNPIDs.txt.txt
