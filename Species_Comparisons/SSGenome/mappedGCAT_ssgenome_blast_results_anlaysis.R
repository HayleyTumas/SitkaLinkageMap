#Find matches between mapped white spruce gene catalog sequences and Sitka spruce geneome
#column names for all results
colnames <- c("qseqid", "sseqid", "sacc","qlen", "slen","pident", "nident", "length", "mismatch", "gaps", "bitscore", "evalue")

#snp chip info
chipkey <- read.csv(paste(datafiles, "SNP_Chip.csv", sep=""))
radkey <- read.csv(paste(datafiles, "RAD_GCAT_firstlook_e11p95.csv", sep=""))
map <- read.csv(paste(datafiles, "RADSNP_SCLM2_Cook4NARC_Map_SNPID.csv", sep=""))


#Load BLAST output
#forward
gcat <- read.csv(paste(datafiles, "GCAT_SSGenome.csv", sep=""), header=F)
#reverse
revgcat <- read.csv(paste(datafiles, "GCAT_SSGenome_REV.csv", sep=""), header=F)
colnames(gcat) <- c(colnames)
colnames(revgcat) <- c(paste("Rev",colnames, sep="_"))

#STEP 1: Find matches that meet thresholds based on e-value and % identity
#Combine forward and reverse results
#create unique match ID
gcat$MatchID <- paste(gcat$qseqid, gcat$sseqid, sep=".")
revgcat$MatchID <- paste(revgcat$Rev_sseqid, revgcat$Rev_qseqid, sep=".")
gcatall <- merge(gcat, revgcat, by="MatchID")#313101
#Filter data 
gcat95e100 <- gcatall[which(gcatall$pident>95 & gcatall$Rev_pident>95 & gcatall$evalue<1e-100 & gcatall$Rev_evalue<1e-100),]#24133
length(unique(gcat95e100$qseqid))
#5237 - Duplicated matches
#Remove duplicated matches
gcat95e100nodup <- gcat95e100[-which(duplicated(gcat95e100$MatchID)),]
#12860
#Remove duplicated GCAT sequences
#i.e. any GCAT sequences that matched to multiple genome scaffolds
#Would be difficult to tell, after filtering for e-value and pident, which match is best
#so just remove any match with a duplicated GCAT sequence
nodups <- gcat95e100nodup[-which(gcat95e100nodup$qseqid %in% gcat95e100nodup[duplicated(gcat95e100nodup$qseqid),]$qseqid),]
write.csv(nodups, file=paste(output, "mappedGCAT_SSGenome_e100p95.csv", sep=""))
#Find any genome scaffolds with multiple GCAT matches
#This would mean that multiple mapped SNPS are on this scaffold and these can be used for validation
dupgenome <- nodups[which(nodups$sseqid %in% nodups[duplicated(nodups$sseqid),]$sseqid),]
#65 duplicated scaaffolds
write.csv(dupgenome, file=paste(output, "Duplicated_Genome_Scaffolds.csv", sep=""))

##STEP 2: Look at map locations for SNPs on duplicated scaffolds
mapchipkey <- chipkey[which(chipkey$SNP %in% map$SNPID),]
snpgen <- merge(dupgenome[,2:3], mapchipkey[,c("SNP", "GCAT")], by.x="qseqid", by.y="GCAT")
radgen <- merge(dupgenome[,2:3], radkey[,c("SNPID", "GCAT")], by.x="qseqid", by.y="GCAT")
colnames(snpgen)[3] <- "SNPID"
snpradgen <- rbind(snpgen, radgen)#65 unique GCAT in 70 rows - some GCAT repeated
mappedgen <- merge(map, snpradgen, by="SNPID")#
write.csv(mappedgen, file=paste(output, "Mapped_SSGenome_Dups.csv", sep=""))
#Link to map positions
#Add GCAT linked to RAD and SNP to map dataframe
key <- data.frame("SNPID"=c(chipkey$SNP, radkey$SNPID), "GCAT"=c(chipkey$GCAT, radkey$GCAT))
mapgcat <- merge(key, map, by="SNPID")#5758-4595 mapped chip SNPS, 1163 RAD with GCAT match
#now get genome scaffold
genmap <- merge(mapgcat, dupgenome[,c("qseqid", "sseqid", "pident", "evalue")], by.x="GCAT", by.y="qseqid")
write.csv(genmap, file=paste(output, "Mapped_ChipRAD_SSGenome_Dups.csv", sep=""))

#NOTE: Examined mismatches by hand due to few duplicates (32)
