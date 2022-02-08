#Prepping data for LepMap3
#set path to DataFiles folder 
datafiles <- #path/

##STEP 1: Get a single SNP per locus.

#up to 3 unique SNPs per locus
#want to get one SNP per locus using SNP with the most reads (determined from vcfTools --site-depth output)
#load site depth output 
depthtab <- read.table(paste(datafiles, "Comb2_SNPS.ldepth", sep=""), sep='\t', header=T)
#create csv file so that is it easier to work with
write.csv(depthtab, file=paste(datafiles, "Comb2_locus_depth.csv", sep=""))
#read back in
depth <- read.csv(paste(datafiles, "Comb2_locus_depth.csv", sep=""), header=T)
#Need to get SNPID associated with the CHR and POS of RADSeq SNPs
snpinfo <- read.csv(paste(datafiles, "RADSEQ_snps.csv", sep=""))
snpinfo$ID <- paste(snpinfo$V1, snpinfo$V2, sep="_")
depth$ID <- paste(depth$CHROM, depth$POS, sep="_")
merged <- merge(snpinfo, depth, by="ID", all.x=FALSE)
ordered <- merged[order(-merged$SUM_DEPTH),]
remdup <- ordered[-which(duplicated(ordered$POS)),]
#get unique SNPs
bestsnps <- remdup[,c("V3", "CHROM", "POS", "SUM_DEPTH", "SUMSQ_DEPTH")]
colnames(bestsnps)[1] <- "SNPID"
write.csv(bestsnps, file=paste(datafiles, "Comb2_BEST_SINGLE_SNP.csv"))
bestsnpid <- as.data.frame(bestsnps$SNPID)
write.table(bestsnpid, file=paste(datafiles, "BEST_COMB2_SNP.txt", sep=""), sep = "\t", col.names = F, row.names = F, quote = F)


##STEP 2: Create pedigree file for LepMap3

#load package
library(vcfR)

#Load merged family vcf from bcftools
bcf <- read.vcfR(file=paste(datafiles, "RAD_comb2_QC.vcf", sep=""))
bcf_inds <- as.data.frame(colnames(bcf@gt))
write.csv(bcf_inds, file=paste(datafiles, "BCF_vcf_Inds.csv", sep=""))

#transpose pedigree file for LepMap3
#load "vertical" pedigree
vped <- read.csv(paste(datafiles, "BCF_vcf_Inds.csv", sep=""), header=F)

ped_t <- as.data.frame(t(vped))

tped_vars <- data.frame(matrix(ncol=2,nrow=6))
tped_vars[,1] <- rep("CHR", 6)
tped_vars[,2] <- rep("POS",6)

tped <- as.data.frame(cbind(tped_vars, ped_t))

write.table(tped, file=paste(datafiles, "BCF_pedtR.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

