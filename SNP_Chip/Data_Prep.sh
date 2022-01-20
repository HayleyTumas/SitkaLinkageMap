
#Note: these files have CHR and POS the same for each SNP
#see R code: all SNPs were assigned CHR and POS using SNP file from GS samples
#CHR is 1 across all SNPs and POS is just the row number

#set environment to run PLINK
conda activate PLINK
#this environment has Plink in it

cd #DataFiles path

#STEP 1: Convert individual family GenomeStudio files to vcf files
#plink files were prepped from GenomeStudio files in R (illumina_to_lgen.R)
#using this tutorial:http://zzz.bwh.harvard.edu/plink/tutorial.shtml
plink --lfile LM1_Chrkey --recode --out LM1P

plink --file LM1P --recode vcf --out LM1_NoZeroSNP

plink --lfile LM2_Chrkey --recode --out LM2P

plink --file LM2P --recode vcf --out LM2_NoZeroSNP

#STEP 2: Filter based on missing data, missingness, and MAF
#create missingness files in Plink, load into R for histogram/diagnostic, use Plink to filter
#need to get plink files by recoding the vcf I had already made in plink
plink --vcf LM1_NoZeroSNP.vcf --recode --out LM1
#now try
plink --file LM1
#that worked
#now make binary ped file to speed up analyses and save space
plink --file LM1 --make-bed --out LM1
#reload the data in binary format using bfile
plink --bfile LM1
#summary stats on missing data
plink --bfile LM1 --missing --out miss_stat
#outputs 2 files with missingness for individuals (imiss) and for SNPs (lmiss)
#now get minor allele frequencies
plink --bfile LM1 --freq --out freq_stat
#now I will load these files into R and take a look at them
#missing maximum for inds is 0.1, for SNPS is 0.2, and MAF max is 0.2
#change ind max missing to 0.15 to match LM2
plink --bfile LM1 --geno 0.2 --mind 0.15 --maf 0.2 --recode vcf --out LM1QC
#double check
plink --vcf LM1QC.vcf --recode --out LM1QC
plink --file LM1QC --missing --out miss_statQC

#LM2
plink --vcf LM2_NoZeroSNP.vcf --recode --out LM2
#now make binary ped file to speed up analyses and save space
plink --file LM2 --make-bed --out LM2
#summary stats on missing data
plink --bfile LM2 --missing --out miss_statL2
#outputs 2 files with missingness for individuals (imiss) and for SNPs (lmiss)
#now get minor allele frequencies
plink --bfile LM2 --freq --out freq_stat
#now I will load these files into R and take a look at them
#missing maximum for inds is 0.15, for SNPS is 0.2, and MAF max is 0.2
plink --bfile LM2 --geno 0.2 --mind 0.15 --maf 0.2 --recode vcf --out LM2QC
#double check in R
plink --vcf LM2QC.vcf --recode --out LM2QC
plink --file LM2QC --missing --out miss_statQC2
#The numbers match so filtering appears to have worked successfully

#STEP 3: merge data to create file with both families for mapping
plink --vcf LM1QC.vcf --recode --out LM1_QC
plink --vcf LM2QC.vcf --recode --out LM2_QC

plink --file LM1_QC --merge LM2_QC.ped LM2_QC.map --recode vcf --out Comb_QC

#STEP 4: Exclude replicated individuals
conda activate VCFTOOLS
#this environment contains VCFtools

vcftools --vcf Comb_QC.vcf --remove SNP_IND_EXCLUDE.txt --out Comb_QC_noreps --recode
