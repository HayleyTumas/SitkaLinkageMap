#Analyze results from blast of mapped GCAT to mapped limber pine sequences

datafiles <- #Path to DataFiles
output <- #Path to Output

#column names for all results
colnames <- c("qseqid", "sseqid", "sacc","qlen", "slen","pident", "nident", "length", "mismatch", "gaps", "bitscore", "evalue")

#snp chip info
chipkey <- read.csv(paste(datafiles, "SNP_Chip.csv", sep=""))
radkey <- read.csv(paste(datafiles, "RAD_GCAT_firstlook_e11p95.csv", sep=""))
map <- read.csv(paste(datafiles, "RADSNP_SCLM2_Cook4NARC_Map_SNPID.csv", sep=""))
pfmap <- read.csv(paste(datafiles, "pflexilis_map.csv", sep=""))

#Load BLAST results
#forward
gcat <- read.csv(paste(datafiles, "OurMap_LimberPine.csv", sep=""), header=F)#3028
#reverse
revgcat <- read.csv(paste(datafiles, "OurMap_LimberPine_REV.csv", sep=""), header=F)#3031
colnames(gcat) <- c(colnames)
colnames(revgcat) <- c(paste("Rev",colnames, sep="_"))

#STEP 1: Identify matches that pass thresholds based on e-value and % identity
#lets look at distribution
hist(gcat$evalue)
nrow(gcat[which(gcat$evalue<1e-100),])#2535 - full list not that long so this is still a large prop
nrow(gcat[which(gcat$pident>95),])#303 -> ok much smaller. Ah well
nrow(gcat[which(gcat$pident>90),])#1721 -> could switch, but will go with 95 for now
#combine forward and reverse files
#create unique match ID
gcat$MatchID <- paste(gcat$qseqid, gcat$sseqid, sep=".")
revgcat$MatchID <- paste(revgcat$Rev_sseqid, revgcat$Rev_qseqid, sep=".")
gcatall <- merge(gcat, revgcat, by="MatchID")#3650
#Try with same thresholds as Norway Spruce
gcat95e100 <- gcatall[which(gcatall$pident>95 & gcatall$Rev_pident>95 & gcatall$evalue<1e-100 & gcatall$Rev_evalue<1e-100),]
#251
#Remove duplicate matches
gcat95e100nodup <- gcat95e100[-which(duplicated(gcat95e100$MatchID)),]
#248 
length(unique(gcat95e100nodup$qseqid))
#248 - only unique GCAT sequences in this subset 
length(unique(gcat95e100nodup$sseqid))
#247 - one limber pine sequences represented twice so will need to remove
nodups <- gcat95e100nodup[-which(gcat95e100nodup$sseqid %in% gcat95e100nodup[duplicated(gcat95e100nodup$sseqid),]$sseqid),]
#246 matches
write.csv(nodups, file=paste(output, "limber_sitka_p95e100.csv", sep=""))
#Due to low number of matches with more rigorous % identity, drop threshold to 90%
#This would be expected due to spruce-pine comparison rather than spruce-spruce comparison
gcat90e100 <- gcatall[which(gcatall$pident>90 & gcatall$Rev_pident>90 & gcatall$evalue<1e-100 & gcatall$Rev_evalue<1e-100),]
#1675 matches
#Remove duplicate matches (i.e. SEQ3-SEQ4 appearing twice)
gcat90e100nodup <- gcat90e100[-which(duplicated(gcat90e100$MatchID)),]
#1490 matches
length(unique(gcat90e100nodup$qseqid))
#1475 - some duplicated GCAT sequences
length(unique(gcat90e100nodup$sseqid))
#1451 - some duplicatd limber pine sequences
#Remove duplicated sequence matches (i.e. SEQ3-SEQ4 and SEQ3-SEQ5 --any matches containing SEQ3 removed)
#Do not trust high quality matches containing a GCAT or limebr pine sequence that is in mulitple high quality matches
noduppf <- gcat90e100nodup[-which(gcat90e100nodup$sseqid %in% gcat90e100nodup[duplicated(gcat90e100nodup$sseqid),]$sseqid),]
#1413 matches
nodupps <- noduppf[-which(noduppf$qseqid %in% noduppf[duplicated(noduppf$qseqid),]$qseqid),]
#1397 matches
write.csv(nodupps, file=paste(output, "limber_sitka_p90e100.csv", sep=""))

##STEP 2: Compare map positions
#first need to get a map that has the GCAT linked to RAD and SNP
key <- data.frame("SNPID"=c(chipkey$SNP, radkey$SNPID), "GCAT"=c(chipkey$GCAT, radkey$GCAT))
colnames(map)[2:3] <- c("LG.PS", "Position.PS")
mapgcat <- merge(key, map[,c("SNPID", "LG.PS", "Position.PS")], by="SNPID")#5758-4595 mapped chip SNPS, 1163 RAD with GCAT match
#Link to BLAST results
genmap90 <- merge(mapgcat, nodupps[,c("qseqid", "sseqid", "pident", "evalue")], by.x="GCAT", by.y="qseqid")
#1514
#Link limber pine to these results
comp90 <- merge(genmap90, pfmap[,c("Linkage.Group", "Composite.position_g9.510", "Transcript.ID")], by.x="sseqid", by.y="Transcript.ID")
#1514
write.csv(comp90, file=paste(output, "Limber_Sitka_p90e100_compare.csv", sep=""))

comp90 <- read.csv(paste(datafiles, "Limber_Sitka_p90e100_compare.csv", sep=""))


#Visually compare 
library(ggplot2)

comp90$LG.PS <- factor(comp90$LG.PS, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
PS_PF_p90 <- ggplot(comp90, aes(x=Position.PS, y=Composite.position_g9.510))+geom_point(size=1)+
  ylab("Position in Limber Pine")+ xlab("Position in Sitka Spruce")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_grid(Linkage.Group~LG.PS)
