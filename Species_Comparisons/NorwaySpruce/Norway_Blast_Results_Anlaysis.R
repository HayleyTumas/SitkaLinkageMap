datafiles <- #Path to DataFiles
output <- #Path to output


#column names for all results
colnames <- c("qseqid", "sseqid", "sacc","qlen", "slen","pident", "nident", "length", "mismatch", "gaps", "bitscore", "evalue")

#Load map data files
chipkey <- read.csv(paste(datafiles, "SNP_Chip.csv", sep=""))
radkey <- read.csv(paste(datafiles, "RAD_GCAT_firstlook_e11p95.csv", sep=""))
map <- read.csv(paste(datafiles, "RADSNP_SCLM2_Cook4NARC_Map_SNPID.csv", sep=""))
pamap <- read.csv(paste(datafiles, "Bernhardssonetal2019.csv", sep=""))

#Load BLAST output
#forward
gcat <- read.csv(paste(datafiles, "OurMap_NorwaySpruce.csv", sep=""), header=F)
#reverse
revgcat <- read.csv(paste(datafiles, "OurMap_NorwaySpruce_REV.csv", sep=""), header=F)
colnames(gcat) <- c(colnames)
colnames(revgcat) <- c(paste("Rev",colnames, sep="_"))

##STEP 1: Filter for e-value and % identity to find matches
#lets look at data distribution
hist(gcat$evalue)
nrow(gcat[which(gcat$evalue<1e-100),])#6640
nrow(gcat[which(gcat$pident>95),])#8791
#combine forward and reverse
#create unique match ID
gcat$MatchID <- paste(gcat$qseqid, gcat$sseqid, sep=".")
revgcat$MatchID <- paste(revgcat$Rev_sseqid, revgcat$Rev_qseqid, sep=".")
gcatall <- merge(gcat, revgcat, by="MatchID")#66631
gcat95e100 <- gcatall[which(gcatall$pident>95 & gcatall$Rev_pident>95 & gcatall$evalue<1e-100 & gcatall$Rev_evalue<1e-100),]
#8174
#Remove duplicated matches
gcat95e100nodup <- gcat95e100[-which(duplicated(gcat95e100$MatchID)),]
#2271 
length(unique(gcat95e100nodup$qseqid))
#2095
length(unique(gcat95e100nodup$sseqid))
#2187
#Since we are using genome scaffolds, there could be multiple gcats to each scaffold
#So will only eliminate GCATs with multiple matches not Norway spruce scaffolds with multiple matches
nodups <- gcat95e100nodup[-which(gcat95e100nodup$qseqid %in% gcat95e100nodup[duplicated(gcat95e100nodup$qseqid),]$qseqid),]
#1935
length(unique(nodups$sseqid))#1873
write.csv(nodups, "mappedGCAT_mappedNorwaySpruceGenome_e100p95.csv", sep=""))

##STEP 2: Match up map positions based on matches
#first need to get a map that has the GCAT linked to RAD and SNP
key <- data.frame("SNPID"=c(chipkey$SNP, radkey$SNPID), "GCAT"=c(chipkey$GCAT, radkey$GCAT))
colnames(map)[2:3] <- c("LG.PS", "Position.PS")
mapgcat <- merge(key, map[,c("SNPID", "LG.PS", "Position.PS")], by="SNPID")#5758-4595 mapped chip SNPS, 1163 RAD with GCAT match
#now get the Norway spruce genome scaffold
genmap <- merge(mapgcat, nodups[,c("qseqid", "sseqid", "pident", "evalue")], by.x="GCAT", by.y="qseqid")
#2071
#now to link to the Bernhardsson et al 2019 map
comp <- merge(genmap, pamap, by.x="sseqid", by.y="Scaffold")
#3234
write.csv(comp, file=paste(output, "NorwaySpruce_Sitka_map_compare.csv", sep=""))

#visual comparison
library(ggplot2)
library(ggpubr)

comp$LG.PS <- factor(comp$LG.PS, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
PS_PA <- ggplot(comp, aes(x=Position.PS, y=Consensus))+geom_point(size=1)+
  ylab("Position in Norway Spruce")+ xlab("Position in Sitka Spruce")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_grid(LG~LG.PS)

length(unique(comp$sseqid))#1873
length(unique(comp$GCAT))#1935
