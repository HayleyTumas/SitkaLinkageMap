datafiles <- #Path to DataFiles
output <- #Path to Output

#column names for all results
colnames <- c("qseqid", "sseqid", "qlen", "slen","pident", "nident", "length", "mismatch", "gaps", "bitscore", "evalue")

#snp chip info
snps <- read.csv(paste(datafiles, "SNP_Chip.csv", sep=""))

#Load data from the forward and reverse BLAST runs
#forward
gcat <- read.csv(paste(datafiles, "RADSeq_FullGCAT.csv", sep=""), header=F)
#reverse
revgcat <- read.csv(paste(datafiles, "RADSeq_FullGCAT_REV.csv", sep=""), header=F)
colnames(gcat) <- c(colnames)
colnames(revgcat) <- c(paste("Rev",colnames, sep="_"))

#STEP 1: Combine forward and reverse results
#create unique match ID
gcat$MatchID <- paste(gcat$qseqid, gcat$sseqid, sep=".")
revgcat$MatchID <- paste(revgcat$Rev_sseqid, revgcat$Rev_qseqid, sep=".")
gcatall <- merge(gcat, revgcat, by="MatchID")#1352
#NOTE: matches are in BLAST multiple times (i.e. SEQ1-SEQ3 match may be listed multiple times)
#AND sequences may be part of multiple matches (i.e. SEQ1-SEQ3 AND SEQ1-SEQ5)
#this may be true if GCAT appears multiple times - could be multiple SNPs on one sequence
#if RADSeq SNP appears multiple times, this indicates a poor match

##STEP 2: Filter based on percent identity and e-value
#Look at histograms to determine good thresholds for determining a reliable match
#base this on percent identity and e-value
hist(gcatall$pident)
hist(gcatall$Rev_pident)
#Determined e value of 1e-11 and % identity of 95 would be good thresholds
#Must pass these thresholds for both forward and reverse blast result
evalue11 <- gcatall[which(gcatall$evalue<1e-11 & gcatall$Rev_evalue<1e-11),]
e11p95 <- evalue11[which(evalue11$pident>95 & evalue11$Rev_pident>95),]
#Remove any duplicated matches (i.e. SEQ1-SEQ3 and SEQ1-SEQ3 multiples)
e11p95nodup <- e11p95[-which(duplicated(e11p95$MatchID)),]
dups <- e11p95nodup[which(e11p95nodup$qseqid %in% e11p95nodup[duplicated(e11p95nodup$qseqid),]$qseqid),]

#STEP 3: Remove duplicated RADSeq sequences and save output
#Because the RADSeq sequence is so short (45 bases)
#the RADSeq sequence can have 100% matches with multiple GCAT sequences
#these 2 high quality matches only occur in 62 RAD loci 
#So exclude any RADSeq loci with multiple matches and the 3 replicate matches
e11p95nodup <- e11p95[-which(duplicated(e11p95$MatchID)),]
e11p95nomulti <- e11p95nodup[-which(e11p95nodup$qseqid %in% dups$qseqid),]
write.csv(e11p95nomulti,paste(output, "RAD_GCAT_firstlook_e11p95.csv", sep=""))
#NOTE: "First look" used in file name before this was determined to be best thresholds
#Keeping this in name to match any downstream analyses

#STEP 4: Create a fasta file for BLASTs against other species
#By creating a fasta with all mapped GCAT sequences
#Will use GCAT sequences for BLASTS because all SNP Chip SNPs are based on GCAT sequences
#And GCAT sequences are much longer than RADSeq sequences
#So need a fasta with GCAT sequences on SNP Chip + GCAT sequences matches to RADSeq SNPs in above file
#RADSeq GCAT key
radgcat <- read.csv(paste(datafiles, "RAD_GCAT_firstlook_e11p95.csv", sep=""))
#sseqid is in different format of SNP Chip GCAT names
#there is a decimal and number at the end that I will need to remove to make them comparable
radgcat$GCAT <- gsub("\\..*","",radgcat$sseqid)
#SNP Chip GCAT key
snpgcat <- read.csv(paste(datafiles, "SNP_Chip.csv", sep=""))
#Map file to get only those SNP Chip SNPs that appear on the final map
#(RADSeq is already subset to only those mapped by using only mapped SNPs in BLAST)
map <- read.csv(paste(datafiles, "RADSNP_SCLM2_Cook4NARC_Map_SNPID.csv", sep=""))
mapsnpgcat <- snpgcat[which(snpgcat$SNP %in% map$SNPID),]
#Get a list of the unique sequences that were mapped across both SNP sets
allgcat <- c(radgcat$GCAT, mapsnpgcat$GCAT)
unigcat <- unique(allgcat)
#Then subset GCAT fasta file 
library(seqinr)
#read in fasta file 
gcat <- read.fasta(file="GCAT_WS-3.3.cluseq.fa")
#Gene Catalog sequence names differ from names in SNP files
#by a decimal followed by a number
#So need to rename sequences in fasta to get rid of .X in GCAT sequence names
names(gcat) <-gsub("\\..*","",names(gcat))
#Subset by those mapped across the SNP chip and RADSeq datasets
map_gcat <- gcat[c(which(names(gcat) %in% unigcat))]
#write the fasta file
write.fasta(sequences = map_gcat, names = names(map_gcat), file.out = "GCAT_Mapped_RADChip.fasta")

##STEP 5: Look at map positions of duplicated GCAT for validation
#Get map positions of mapped GCAT sequences for both RADSeq and SNP Chip
snpmerge <- merge(map, snpgcat[,c("SNP", "GCAT")], by.x="SNPID", by.y="SNP")
radmerge <- merge(map, radgcat[,c("SNPID", "GCAT")], by="SNPID")
#combine the 2 datasets
gcatmap <- rbind(snpmerge, radmerge) #only 5758
write.csv(gcatmap, file=paste(output, "GCAT_positions_RADChip.csv", sep=""))
#Find GCATs with multiple SNPs
fullsame <- gcatmap[which(gcatmap$GCAT %in% gcatmap[duplicated(gcatmap$GCAT),]$GCAT),]
#670, (326 unique GCATs with 2-3 SNPs)
#Now determine how many of those do not fall on the same linkage group
#separate out by those that have only 2 SNPs and those that have 3 SNPs
fullsamedup1 <- fullsame[-which(duplicated(fullsame$GCAT)),]
fullsamedups <- fullsame[which(duplicated(fullsame$GCAT)),]
fullsamedup2 <- fullsamedups[-which(duplicated(fullsamedups$GCAT)),]
fullsamedup3 <- fullsamedups[which(duplicated(fullsamedups$GCAT)),]
#GCAT with 2 SNPs
match2 <- fullsamedup1[-which(fullsamedup1$GCAT %in% fullsamedup3$GCAT),]
match2dups <- merge(match2, fullsamedup2, by="GCAT", all.y=F)
#How many are on different linkage groups?
mismatch2 <- match2dups[-which(match2dups$LG.x==match2dups$LG.y),]
#26
#GCAT with 3 SNPs
firststep <- merge(fullsamedup1[which(fullsamedup1$GCAT %in% fullsamedup3$GCAT),], fullsamedup2[which(fullsamedup2$GCAT %in% fullsamedup3$GCAT),], by="GCAT")
colnames(fullsamedup3) <- paste(colnames(fullsamedup3), "3", sep=".")
match3dups <- merge(firststep, fullsamedup3, by.x="GCAT", by.y="GCAT.3")
#How many are on different linkage groups?
mismatch3 <- match3dups[-which(match3dups$LG.x==match3dups$LG.y&match3dups$LG.x==match3dups$LG.3),]
#2 - only 3rd SNP does not match
#what is the percentage of GCAT with SNPs on the same linkage groups?
#markers
(670-((26*2)+2))/670
(326-(26+2))/326
#Create files with same and different linkage groups
#create file with SNPs on different linkage groups
first <- merge(fullsamedup1, fullsamedup2, by="GCAT")
second <- merge(first, fullsamedup3, by.x="GCAT", by.y="GCAT.3", all.x=T)
fullmismatches <- second[which(second$GCAT %in% mismatch2$GCAT | second$GCAT %in% mismatch3$GCAT),]
write.csv(fullmismatches, file=paste(output, "Map_Mismatch.csv", sep=""))
#creat file with SNPs on the same linkage group
matches_first <- fullsame[-which(fullsame$GCAT %in% fullmismatches$GCAT),]
matchdup1 <- matches_first[-which(duplicated(matches_first$GCAT)),]
matchdups <- matches_first[which(duplicated(matches_first$GCAT)),]
matchdup2 <- matchdups[-which(duplicated(matchdups$GCAT)),]
matchdup3 <- matchdups[which(duplicated(matchdups$GCAT)),]
matchmerge1 <- merge(matchdup1, matchdup2, by="GCAT")
colnames(matchdup3) <- paste(colnames(matchdup3), "3", sep=".")
matches <- merge(matchmerge1, matchdup3, by.x="GCAT", by.y="GCAT.3", all.x=T)
matches$Diff_Pos_Set1 <- abs(matches$Position.x-matches$Position.y)
matches$Diff_Pos_Set2 <- abs(matches$Position.y-matches$Position.3)
write.csv(matches, file=paste(output, "Map_Matches.csv", sep=""))

#STEP 6: Visualize validation
#visual
library(ggplot2)
library(ggpubr)

match2dups$LG.x<- factor(match2dups$LG.x, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
match2dups$LG.y<- factor(match2dups$LG.y, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))

gcat_dups <- ggplot(match2dups, aes(x=Position.x, y=Position.y))+geom_point(size=1)+
  ylab("Position in First Duplicate")+ xlab("Position in Second Duplicate")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_grid(LG.y~LG.x)

##STEP 7: Calcualte validation correlation
sameLG <- match2dups[-which(match2dups$GCAT %in% mismatch2$GCAT),]
Match2LGs <- split(sameLG , f = sameLG$LG.x )
Match2corr <- list()

for(x in 1:12){
  nam <- paste("LG", x, sep="")
  Match2corr[[nam]]  <- data.frame("Method"=cor.test(Match2LGs[[x]]$Position.x, Match2LGs[[x]]$Position.y, method="kendall", ci=FALSE)$method, "Estimate"=cor.test(Match2LGs[[x]]$Position.x, Match2LGs[[x]]$Position.y, method="kendall", ci=FALSE)$estimate, "Pvalue"=as.numeric(cor.test(Match2LGs[[x]]$Position.x, Match2LGs[[x]]$Position.y, method="kendall", ci=FALSE)$p.value))
}

Match2corrs <- do.call("rbind", Match2corr) #computing p value impossible with ties
mean(abs(Match2corrs$Estimate))
#0.9380623
min(Match2corrs$Estimate)
# 0.7673883
max(Match2corrs$Estimate)
#1
