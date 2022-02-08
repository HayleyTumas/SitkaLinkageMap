#Path to DataFiles folder
datafiles <- #pathtoDataFiles/

##Step 1: Convert from GenomeStudio file to lgen for uptake into Plink

library(dplyr)
library(tidyr)
library(argparse)
require(readxl)
require(tidyverse)
require(tibble)

#FinalReport creation in GenomeStudio:
#Analysis->Reports->Wizard->Final Report
#Next to select all samples groups
#Remove excluded samples from report
#Remove zeroed SNPs, Include intensity only SNPs
#Remove positive control samples from Family 1 (C3)
#Standard format
#Columns: SNP Name, Sample ID, Allele1 - Top, Allele2 - Top, Chr, Position 
#*remove GC Score and add Chr and Position 
#Comma separated file
#Then follow following steps to convert to format to be input into PLINK

#need to create "universal" chromosome/position for SNPs for combining vcfs
GSfile <- read.csv(paste(datafiles, file="GS_FinalReport.csv", sep=""),skip = 9, header = T)
SNPIDs <- unique(GSfile$SNP.Name)
key <- data.frame("SNP.Name"=c(SNPIDs), "Chr"=c(rep(1, length(SNPIDs))), "Position"=c(1:length(SNPIDs)))
write.csv(key, file="SNP_ChrPos_Key_Combine.csv")
#merge it with final report for LM1 and LM2
#LM1 (C3)
Finalreport<-read.csv(file="C3_NoPosCon_FinalReport_NoZeroSNP.csv",skip = 9, header = T)
#Fixing sample ID name
Finalreport <- Finalreport[,-c(5,6)]
Finalreport$Sample.ID <- sub(".*# *(.*?) *#.*", "\\1", Finalreport$Sample.ID)
Finalreportm <- merge(Finalreport, key, by="SNP.Name")
Finalreport <- as.data.frame(Finalreportm)
Finalreport<-filter(Finalreport, Allele1...Top != 'I')
Finalreport["empty"]="0"
Finalreport['fid']<-Finalreport$Sample.ID
map<-select(Finalreport, Chr, SNP.Name, empty, Position)
map<-map[!duplicated(map),]
map<-map[complete.cases(map),]
lgen<-select(Finalreport, fid, Sample.ID, SNP.Name, Allele1...Top, Allele2...Top)
lgen<-lgen[!duplicated(lgen),]
lgen<-lgen[complete.cases(lgen),]
lgen<-filter(lgen, Allele1...Top != "-" & Allele2...Top != "-")
write.table(map, file = paste(datafiles, "LM1_Chrkey.map", sep=""), sep = "\t", col.names = F, row.names = F, quote = F)
write.table(lgen, file =paste(datafiles, "LM1_Chrkey.lgen", sep=""), sep = "\t", col.names = F, row.names = F, quote = F)

#make fam file
Finalreport["Father"] <- "0"
Finalreport["Mother"] <- "0"
Finalreport["Sex"] <- "0"
Finalreport["Phenotype"] <- "0"
fam <- select(Finalreport, fid, Sample.ID, Father, Mother, Sex, Phenotype)
fam <- fam[!duplicated(fam$Sample.ID),]
write.table(fam, file = paste(datafiles, "LM1_Chrkey.fam", sep=""), sep = "\t", col.names = F, row.names = F, quote = F)

#LM2 (C1)
#clean environment and reload key
key <- read.csv(paste(datafiles, "SNP_ChrPos_Key_Combine.csv", sep=""))
#using file without zeroed SNPs
Finalreport<-read.csv(paste(datafiles, file="C1_FinalReport_NoZeroSNPs_Standard.csv", sep=""),skip = 9, header = T)
#Fixing sample ID name
Finalreport <- Finalreport[,-c(5,6)]
Finalreport$Sample.ID <- sub(".*# *(.*?) *#.*", "\\1", Finalreport$Sample.ID)
Finalreportm <- merge(Finalreport, key, by="SNP.Name")
Finalreport <- as.data.frame(Finalreportm)
Finalreport<-filter(Finalreport, Allele1...Top != 'I')
Finalreport["empty"]="0"
Finalreport['fid']<-Finalreport$Sample.ID
map<-select(Finalreport, Chr, SNP.Name, empty, Position)
map<-map[!duplicated(map),]
map<-map[complete.cases(map),]
lgen<-select(Finalreport, fid, Sample.ID, SNP.Name, Allele1...Top, Allele2...Top)
lgen<-lgen[!duplicated(lgen),]
lgen<-lgen[complete.cases(lgen),]
lgen<-filter(lgen, Allele1...Top != "-" & Allele2...Top != "-")
write.table(map, file = paste(datafiles, "LM2_Chrkey.map", sep=""), sep = "\t", col.names = F, row.names = F, quote = F)
write.table(lgen, file = paste(datafiles, "LM2_Chrkey.lgen", sep=""), sep = "\t", col.names = F, row.names = F, quote = F)

#make fam file
Finalreport["Father"] <- "0"
Finalreport["Mother"] <- "0"
Finalreport["Sex"] <- "0"
Finalreport["Phenotype"] <- "0"
fam <- select(Finalreport, fid, Sample.ID, Father, Mother, Sex, Phenotype)
fam <- fam[!duplicated(fam$Sample.ID),]
write.table(fam, file = paste(datafiles, "LM2_Chrkey.fam", sep=""), sep = "\t", col.names = F, row.names = F, quote = F)



##STEP 2: Determine QC thresholds - see accompanying Plink code

#LM1
#load data 
#individual missing data
ind_miss <- read.table(paste(datafiles, "miss_stat.imiss", sep=""), header=T)
snp_miss <- read.table(paste(datafiles, "miss_stat.lmiss", sep=""), header=T)
maf <- read.table(paste(datafiles, "freq_stat.frq", sep=""), header=T)

hist(ind_miss$F_MISS)
hist(ind_miss$N_MISS)
nrow(ind_miss[which(ind_miss$F_MISS<0.10),])
#314
nrow(ind_miss[which(ind_miss$F_MISS<0.05),])
#304
nrow(ind_miss[which(ind_miss$F_MISS<0.08),])
#312
nrow(ind_miss[which(ind_miss$F_MISS<0.15),])
#315

#go with 0.10 missing data for individuals

hist(snp_miss$F_MISS)
nrow(snp_miss[which(snp_miss$F_MISS<0.25),])
#10861
nrow(snp_miss[which(snp_miss$F_MISS<0.2),])
#10798
nrow(snp_miss[which(snp_miss$F_MISS<0.1),])
#10633
#go with 0.2 missing data for SNPs (genotypes in PLINK)

hist(maf$MAF)
nrow(maf[which(maf$MAF>0.2),])
#use 0.2 since everything else is noninformative for linkage mapping

#filtered double check (filter in plink then recreate stat files)
ind_missQC <- read.table(paste(datafiles, "miss_statQC.imiss", sep=""), header=T)
snp_missQC <- read.table(paste(datafiles, "miss_statQC.lmiss", sep=""), header=T)
maf <- read.table(paste(datafiles, "freq_stat.frq", sep=""), header=T)

hist(ind_miss$F_MISS)

#LM2
#load data 
#individual missing data
ind_miss <- read.table(paste(datafiles, "miss_stat.imiss", sep=""), header=T)
snp_miss <- read.table(paste(datafiles, "miss_stat.lmiss", sep=""), header=T)
maf <- read.table(paste(datafiles, "freq_stat.frq", sep=""), header=T)

hist(ind_miss$F_MISS)
hist(ind_miss$N_MISS)
nrow(ind_miss[which(ind_miss$F_MISS<0.10),])
#296
nrow(ind_miss[which(ind_miss$F_MISS<0.15),])
#305


#go with 0.15 missing data

hist(snp_miss$F_MISS)
nrow(snp_miss[which(snp_miss$F_MISS<0.25),])
#11569
nrow(snp_miss[which(snp_miss$F_MISS<0.2),])
#11460

#go with 0.2 missing data

hist(maf$MAF)
nrow(maf[which(maf$MAF>0.2),])
#4097
#use 0.2 since everything else is noninformative for linkage mapping

#filtered double check (filter in plink then recreate stat files)
ind_missQC <- read.table("miss_statQC.imiss", header=T)
snp_missQC <- read.table("miss_statQC.lmiss", header=T)
hist(ind_miss$F_MISS)


##STEP 3: Get sample IDs from family vcf files. 

library(vcfR)
#get individuals from LM1
lm1 <- read.vcfR(paste(datafiles, "LM1QC.vcf", sep=""))
strwrap(lm1@meta)
lm1@gt[1:10, 1:20]
samples <- as.data.frame(colnames(lm1@gt))
write.csv(samples, file=paste(datafiles, "LM1QC_Inds.csv", sep=""))

#get individuals from LM2
lm2 <- read.vcfR(paste(datafiles, "LM2QC.vcf", sep=""))
strwrap(lm1@meta)
lm1@gt[1:10, 1:20]
samples2 <- as.data.frame(colnames(lm2@gt))

write.csv(samples2, file=paste(datafiles,"LM2QC_Inds.csv", sep=""))

#get all individuals in combined dataset
allsamples <- rbind(samples, samples2)

write.csv(allsamples, file=paste(datafiles,"SNP_2FAMCombine_Inds.csv", sep=""))

##STEP 4: Get list of duplicate individuals.
snp <- read.csv(paste(datafiles,"SNP_2FAMCombine_Inds.csv", sep=""))

snp_reps <- snp[which(duplicated(snp$IND)),]
snp_repind <- snp[which(snp$IND %in% snp_reps$IND),]
write.csv(snp_repind, file=paste(datafiles,"SNP_IND_REPS.csv", sep=""))

#create text file to exclude these reps (the ones withtout the r) from vcf files
exclude_inds <- as.data.frame(c(paste(snp_reps$SNP_Pre, snp_reps$SNP_IND, sep="_")))
write.table(exclude_inds, file=paste(datafiles,"SNP_IND_EXCLUDE.txt", sep=""), sep = "\t", col.names = F, row.names = F, quote = F)

##Step 5: Create pedigree file

#use above single family sample lists to create extra columns for LepMap3
LM1_inds <- data.frame("ID"=c(as.character(samples[-1,])))
#Due to nature of parent names, both parent IDs were listed last
#Reversing this order puts the parents at top of dataframe in rows 1 and 2
LM1_inds$order <- c(1:nrow(LM1_inds))
LM1samples <- LM1_inds[order(-LM1_inds$order),]
#More information on pedigree file found on LepMap3 wiki
#All individuals here are in the same family 
#Parent IDs have 0 for parent ID
#These species are monoecious so does not matter of Parent in row 1 or row 2 are listed as P1 or P2
#Sex randomly assigned to 2 parents and listed as unknown for all offspring
LM1_pedigree <- data.frame("FAM"=c(rep("LM1", nrow(LM1samples))), "IND"=c(LM1samples$ID), "P1"=c(0,0, rep(LM1samples[1,1], nrow(LM1samples)-2)), "P2"=c(0, 0, rep(LM1samples[2,1], nrow(LM1samples)-2)), "Sex"=c(1,2, rep(0, nrow(LM1samples)-2)), "Fill"=c(rep(1, nrow(LM1samples))))

LM2_inds <- data.frame("ID"=c(as.character(samples2[-1,])))
LM2_inds$order <- c(1:nrow(LM2_inds))
LM2samples <- LM2_inds[order(-LM2_inds$order),]
LM2_pedigree <- data.frame("FAM"=c(rep("LM2", nrow(LM2samples))), "IND"=c(LM2samples$ID), "P1"=c(0,0, rep(LM2samples[1,1], nrow(LM2samples)-2)), "P2"=c(0, 0, rep(LM2samples[2,1], nrow(LM2samples)-2)), "Sex"=c(1,2, rep(0, nrow(LM2samples)-2)), "Fill"=c(rep(1, nrow(LM2samples))))

#combine both families for using the combined vcf file
combined <- rbind(LM1_pedigree, LM2_pedigree)
#save ALL individuals before the replicates had been removed
write.table(combined, file=paste(datafiles,"Combine_inds.csv", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep=",")

#Remove replicates from pedigree file and save
ped <- read.csv(paste(datafiles,"Combine_Inds.csv", sep=""), header=F) 
ped_noreps <- ped[-which(ped$V2 %in% exclude_inds$`c(paste(snp_reps$SNP_Pre, snp_reps$SNP_IND, sep = "_"))`),]
write.csv(ped_noreps, file=paste(datafiles,"Combine_Inds_Noreps.csv", sep=""))

#reload file of all individuals
vped <- read.csv(paste(datafiles,"Combine_Inds_Noreps.csv", sep=""), header=F)

#transpose file and add columns for LepMap3 (see wiki)
ped_t <- as.data.frame(t(vped))

tped_vars <- data.frame(matrix(ncol=2,nrow=6))
tped_vars[,1] <- rep("CHR", 6)
tped_vars[,2] <- rep("POS",6)

tped <- as.data.frame(cbind(tped_vars, ped_t))

#save final pedigree file that will be input into LepMap3
write.table(tped, file=paste(datafiles,"Comb_ped_NoReps.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

