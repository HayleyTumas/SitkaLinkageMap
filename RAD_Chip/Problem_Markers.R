datafiles <= #Path to DataFiles/

##Step 1: Identifying markers causing gaps
#Using map from combined family data using LM2 during separate chromosomes
map <- read.csv(paste(datafiles, "RADSNP_SCLM2_Map_SNPID.csv", sep=""), header=T)

#First looked to see on which chromosomes gaps fell
#visualize
lg.levels <- c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12")
map$LG <- factor(map$LG, levels=lg.levels)
lowgapplot <- ggplot(lowgaps, aes(x=Position, y=LG)) + geom_point()
regplot <- ggplot(map, aes(x=Position, y=LG)) + geom_point()

#then determined gap between markers
map_order <- map[order(map$LG, map$Position),]

require(dplyr)
gaps <- map_order %>%
  group_by(LG) %>%
  mutate(Gaps = Position - lag(Position, default = Position[1]))

write.csv(gaps, paste(datafiles, file="RADSNP_SCLM2_Gaps.csv", sep=""))

#then manually trimmed gaps in 1, 2, 4, and 9
#identified which linkage groups based on visual graph
#determined which markers were causing large gap using gap calculate above
#trimmed that marker and any preceding (top of linkage group) or later (end of linkage group) markers
trim <- read.csv(paste(datafiles, "RADSNP_SCLM2_Gaps_Trim.csv", sep=""))
#find trimmed markers
gap_markers <- map[-which(map$SNPID %in% trim$SNPID),]
write.csv(gap_markers, paste(datafiles, file="RADSNP_SCLM2_GapMarkers.csv", sep=""))


##STEP 2: Find any markers that consistently do not group on the same linkage group
#when comparing the map to: family maps (2), chip map, RAD map, and GCAT
#load maps
comb2 <- read.csv(paste(datafiles, "RADSNP_SCLM2_Map_SNPID.csv", sep=""), header=T)
lm1 <- read.csv(paste(datafiles, "RADSNP_LM1_Map_SNPID.csv", sep=""), header=T)
lm2 <- read.csv(paste(datafiles, "RADSNP_LM2_Map_SNPID.csv", sep=""), header=T)
chip <- read.csv(paste(datafiles, "CombQCNR_LM2SC_Map_SNPID.csv", sep=""))
RAD <- read.csv(paste(datafiles, "RADQC_2FAMSSSCLM2_Map_SNPID.csv", sep=""))
ws <- read.csv(paste(datafiles, "SNP_Chip.csv", sep=""), header=T)

#Rename columns for merging
colnames(comb2)[2:5] <- paste(colnames(comb2)[2:5], "Comb2", sep=".")
colnames(lm1)[2:5] <- paste(colnames(lm1)[2:5], "LM1", sep=".")
colnames(lm2)[2:5] <- paste(colnames(lm2)[2:5], "LM2", sep=".")
colnames(chip)[2:5] <- paste(colnames(chip)[2:5], "Chip", sep=".")
colnames(RAD)[2:5] <- paste(colnames(RAD)[2:5], "RAD", sep=".")

#RADSNP vs RADSNP LM1
#Merge the two maps based on SNPID
lm1comb <- merge(comb2, lm1, by="SNPID")
#Then aggregate based on linkage group name to determine which linkage groups align
#i.e. if 200 of the same markers fall on LG1 in Map 1 and LG4 on map 2, then LG1 and LG4 are the same chromosome
lm1comb_lg <- aggregate(lm1comb$Position.LM1, by = lm1comb[c('LG.Comb2', 'LG.LM1')], length)
write.csv(lm1comb_lg, paste(datafiles, file="LM1_Comb2_LGs.csv", sep=""))
#visual look to determine which linkage groups match up
complm1comb_lg <- aggregate(complm1comb$SNPID, by = complm1comb[c('LG.Comb2', 'LG.LM1')], length)
write.csv(complm1comb,paste(datafiles, complm1comb, "LM1_Comb2_Compare.csv", sep=""))
complm1miss <- lm1comb[-which((lm1comb$LG.Comb2=="LG1"&lm1comb$LG.LM1=="LG1")|(lm1comb$LG.Comb2=="LG10"&lm1comb$LG.LM1=="LG9")|(lm1comb$LG.Comb2=="LG11"&lm1comb$LG.LM1=="LG12")|(lm1comb$LG.Comb2=="LG12"&lm1comb$LG.LM1=="LG11")|(lm1comb$LG.Comb2=="LG2"&lm1comb$LG.LM1=="LG2")|(lm1comb$LG.Comb2=="LG3"&lm1comb$LG.LM1=="LG8")|(lm1comb$LG.Comb2=="LG4"&lm1comb$LG.LM1=="LG3")|(lm1comb$LG.Comb2=="LG5"&lm1comb$LG.LM1=="LG4")|(lm1comb$LG.Comb2=="LG6"&lm1comb$LG.LM1=="LG7")|(lm1comb$LG.Comb2=="LG7"&lm1comb$LG.LM1=="LG10")|(lm1comb$LG.Comb2=="LG8"&lm1comb$LG.LM1=="LG6")|(lm1comb$LG.Comb2=="LG9"&lm1comb$LG.LM1=="LG5")),]
write.csv(complm1miss,paste(datafiles,  "LM1_Comb2_LGMiss.csv", sep=""))

#RADSNP vs RADSNP LM2
lm2comb <- merge(comb2, lm2, by="SNPID")
lm2comb_lg <- aggregate(lm2comb$Position.LM2, by = lm2comb[c('LG.Comb2', 'LG.LM2')], length)
write.csv(lm2comb_lg,paste(datafiles,  file="LM2_Comb2_LGs.csv", sep=""))
complm2comb <- lm2comb[-which((lm2comb$LG.Comb2=="LG5"&lm2comb$LG.LM2=="LG2")),]
write.csv(complm2comb,paste(datafiles,  file="LM2_Comb2_Compare.csv", sep=""))
complm2miss <- lm2comb[which((lm2comb$LG.Comb2=="LG5"&lm2comb$LG.LM2=="LG2")),]
write.csv(complm2miss, paste(datafiles, file="LM2_Comb2_LGMiss.csv", sep=""))

#RADSNP vs SNP Chip combined map using LM2 during Separate Chromosomes
radcchipm <- merge(comb2, chip, by="SNPID", all.x=F)
radc_chip_lg <- aggregate(radcchipm$SNPID, by=radcchipm[c('LG.Comb2', 'LG.Chip')], length)
write.csv(radc_chip_lg, paste(datafiles, file="Comb2_ChipCombSCLM2.csv", sep=""))
compComb2Chip <- radcchipm[-which((radcchipm$LG.Comb2=="LG5"&radcchipm$LG.Chip=="LG5")|(radcchipm$LG.Comb2=="LG7"&radcchipm$LG.Chip=="LG5")),]
write.csv(compComb2Chip, paste(datafiles, file="Comb2_Chip_Compare.csv", sep=""))
compComb2Chipmiss <- radcchipm[which((radcchipm$LG.Comb2=="LG5"&radcchipm$LG.Chip=="LG5")|(radcchipm$LG.Comb2=="LG7"&radcchipm$LG.Chip=="LG5")),]
write.csv(compComb2Chipmiss, paste(datafiles, file="Comb2_Chip_LGMiss.csv", sep=""))

#RADSNP vs RADSeq combined map using LM2 during Separate Chromosomes
RADC2 <- merge(comb2, RAD, by="SNPID")
RADC2_lg <- aggregate(RADC2$SNPID, by=RADC2[,c("LG.Comb2", "LG.RAD")], length)
write.csv(RADC2_lg, paste(datafiles, "Comb2_RADSS_SCL2_LG.csv", sep=""))
compradcomb <- RADC2[-which((RADC2$LG.Comb2=="LG4"&RADC2$LG.RAD=="LG11")|(RADC2$LG.Comb2=="LG2"&RADC2$LG.RAD=="LG5")|(RADC2$LG.Comb2=="LG2"&RADC2$LG.RAD=="LG7")|(RADC2$LG.Comb2=="LG4"&RADC2$LG.RAD=="LG7")|(RADC2$LG.Comb2=="LG5"&RADC2$LG.RAD=="LG7")|(RADC2$LG.Comb2=="LG10"&RADC2$LG.RAD=="LG9")),]
write.csv(compradcomb, paste(datafiles, "RADSCLM2_Comb2_Compare.csv", sep=""))
compradmiss <- RADC2[which((RADC2$LG.Comb2=="LG4"&RADC2$LG.RAD=="LG11")|(RADC2$LG.Comb2=="LG2"&RADC2$LG.RAD=="LG5")|(RADC2$LG.Comb2=="LG2"&RADC2$LG.RAD=="LG7")|(RADC2$LG.Comb2=="LG4"&RADC2$LG.RAD=="LG7")|(RADC2$LG.Comb2=="LG5"&RADC2$LG.RAD=="LG7")|(RADC2$LG.Comb2=="LG10"&RADC2$LG.RAD=="LG9")),]
write.csv(compradmiss, paste(datafiles, "RADSCLM2_Comb2_LGMiss.csv", sep=""))

#RADSNP vs White Spruce Map
comb2ws <- merge(comb2, ws, by.x="SNPID", by.y="SNP")
mapped2 <- comb2ws[which(comb2ws$Mapped=="Y"),]
comb2ws_lg <- aggregate(mapped2$SNPID, by=mapped2[c('LG.Comb2', 'Composite.map.linkage.group')], length)
write.csv(comb2ws_lg, paste(datafiles, file="Comb2_WS.csv", sep=""))
compcomb2ws <- mapped2[which((mapped2$LG.Comb2=="LG1"&mapped2$Composite.map.linkage.group=="LG07")|(mapped2$LG.Comb2=="LG10"&mapped2$Composite.map.linkage.group=="LG12")|(mapped2$LG.Comb2=="LG11"&mapped2$Composite.map.linkage.group=="LG10")|(mapped2$LG.Comb2=="LG12"&mapped2$Composite.map.linkage.group=="LG09")|(mapped2$LG.Comb2=="LG2"&mapped2$Composite.map.linkage.group=="LG01")|(mapped2$LG.Comb2=="LG3"&mapped2$Composite.map.linkage.group=="LG08")|(mapped2$LG.Comb2=="LG4"&mapped2$Composite.map.linkage.group=="LG02")|(mapped2$LG.Comb2=="LG5"&mapped2$Composite.map.linkage.group=="LG05")|(mapped2$LG.Comb2=="LG6"&mapped2$Composite.map.linkage.group=="LG03")|(mapped2$LG.Comb2=="LG7"&mapped2$Composite.map.linkage.group=="LG06")|(mapped2$LG.Comb2=="LG8"&mapped2$Composite.map.linkage.group=="LG11")|(mapped2$LG.Comb2=="LG9"&mapped2$Composite.map.linkage.group=="LG04")),]
write.csv(compcomb2ws, paste(datafiles, file="Compare_Comb2_WS.csv", sep=""))
compcomb2wsMiss <- mapped2[-which((mapped2$LG.Comb2=="LG1"&mapped2$Composite.map.linkage.group=="LG07")|(mapped2$LG.Comb2=="LG10"&mapped2$Composite.map.linkage.group=="LG12")|(mapped2$LG.Comb2=="LG11"&mapped2$Composite.map.linkage.group=="LG10")|(mapped2$LG.Comb2=="LG12"&mapped2$Composite.map.linkage.group=="LG09")|(mapped2$LG.Comb2=="LG2"&mapped2$Composite.map.linkage.group=="LG01")|(mapped2$LG.Comb2=="LG3"&mapped2$Composite.map.linkage.group=="LG08")|(mapped2$LG.Comb2=="LG4"&mapped2$Composite.map.linkage.group=="LG02")|(mapped2$LG.Comb2=="LG5"&mapped2$Composite.map.linkage.group=="LG05")|(mapped2$LG.Comb2=="LG6"&mapped2$Composite.map.linkage.group=="LG03")|(mapped2$LG.Comb2=="LG7"&mapped2$Composite.map.linkage.group=="LG06")|(mapped2$LG.Comb2=="LG8"&mapped2$Composite.map.linkage.group=="LG11")|(mapped2$LG.Comb2=="LG9"&mapped2$Composite.map.linkage.group=="LG04")),]
write.csv(compcomb2wsMiss, paste(datafiles, file="Compare_Comb2_WS_LGMiss.csv", sep=""))

#Find SNPs that are causing repeated miss-matches in LGs
CombLM1 <- read.csv(paste(datafiles, "LM1_Comb2_LGMiss.csv", sep=""))
CombLM2 <- read.csv(paste(datafiles, "LM2_Comb2_LGMiss.csv", sep=""))          
CombChip <- read.csv(paste(datafiles, "Comb2_Chip_LGMiss.csv", sep=""))
CombWS <- read.csv(paste(datafiles, "Compare_Comb2_WS_LGMiss.csv", sep=""))
CombRAD <- read.csv(paste(datafiles, "RADSCLM2_Comb2_LGMiss.csv", sep=""))


misses <- data.frame("SNPID"=c(CombLM1[,c("SNPID")], CombLM2[,c("SNPID")], CombChip[,c("SNPID")], CombWS[,c("SNPID")], CombRAD[,c("SNPID")]))
misses$Match <- c(rep("LM1", nrow(CombLM1)), rep("LM2", nrow(CombLM2)), rep("Chip", nrow(CombChip)), rep("WS", nrow(CombWS)), rep("RAD", nrow(CombRAD)))
missnum <- aggregate(misses$Match, by=list(misses$SNPID), length)
max(missnum$x)
#3 - so none are missed across all 5, which they could not be
bigmiss <- missnum[which(missnum$x>2),]
#1
smallmiss <- missnum[which(missnum$x>1),]
#20 -> ok so no RAD are miss grouping across all 3 sets. Guessing its LM1 and RADSCLM2 causing the new 12
write.csv(smallmiss, paste(datafiles, file="MissLG_multi.csv", sep=""))

##STEP 3: Find markers that are causing issues with position
#Did this using Cook's distance of greater than 4/N where N is sample size
#based this on idea that, when comparing 2 maps, positions that are the same or have a high correlation would fall on linear line
#this line could be positive or negative depending on if the linkage group of one map is inverted compared to the other
#So Cook's distance is typically used to identify outliers in linear regression
#so basically did the same thing here.
#used files that only had markers that fell on same linkage group (see above)
#load data
RADC2 <- read.csv(paste(datafiles, "RADSCLM2_Comb2_Compare.csv", sep=""))
LM1C2 <- read.csv(paste(datafiles, "LM1_Comb2_Compare.csv", sep=""))
LM2C2 <- read.csv(paste(datafiles, "LM2_Comb2_Compare.csv", sep=""))
ChipC2 <- read.csv(paste(datafiles, "Comb2_Chip_Compare.csv", sep=""))
WSC2 <- read.csv(paste(datafiles, "Compare_Comb2_WS.csv", sep=""))

#RADSNP vs LM1
LM1LGs <- split( LM1C2 , f = LM1C2$LG.Comb2 )

for(x in 1:12){
  LM1LGs[[x]]$Cook  <- (cooks.distance(lm(Position.LM1~Position.Comb2, LM1LGs[[x]])))
}

LM1Outs <- list()
for(i in 1:12){
  LM1Outs[[i]] <- as.data.frame(subset(LM1LGs[[i]], Cook>(4/length(LM1LGs[[i]]$Cook))))
}
#get all problematic markers
LM1PM <- do.call("rbind", LM1Outs)
#235

#RADSNP vs LM2
LM2LGs <- split( LM2C2 , f = LM2C2$LG.Comb2 )

for(x in 1:12){
  LM2LGs[[x]]$Cook  <- (cooks.distance(lm(Position.LM2~Position.Comb2, LM2LGs[[x]])))
}

LM2Outs <- list()
for(i in 1:12){
  LM2Outs[[i]] <- as.data.frame(subset(LM2LGs[[i]], Cook>(4/length(LM2LGs[[i]]$Cook))))
}
#get all problematic markers
LM2PM <- do.call("rbind", LM2Outs)
#162

#RADSNP vs Chip
CLGs <- split( ChipC2 , f = ChipC2$LG.Comb2 )

for(x in 1:12){
  CLGs[[x]]$Cook  <- (cooks.distance(lm(Position.Chip~Position.Comb2, CLGs[[x]])))
}

COuts <- list()
for(i in 1:12){
  COuts[[i]] <- as.data.frame(subset(CLGs[[i]], Cook>(4/length(CLGs[[i]]$Cook))))
}
#get all problematic markers
CPM <- do.call("rbind", COuts)
#247

#RADSNP vs RADSeq
LGs <- split( RADC2 , f = RADC2$LG.Comb2 )

for(x in 1:12){
  LGs[[x]]$Cook  <- (cooks.distance(lm(Position.RAD~Position.Comb2, LGs[[x]])))
}

Outs <- list()
for(i in 1:12){
  Outs[[i]] <- as.data.frame(subset(LGs[[i]], Cook>(4/length(LGs[[i]]$Cook))))
}
#get all problematic markers
RADPM <- do.call("rbind", Outs)
#417

#RADSNP vs WS
WSLGs <- split( WSC2 , f = WSC2$LG.Comb2 )

for(x in 1:12){
  WSLGs[[x]]$Cook  <- (cooks.distance(lm(Composite.map.position~Position.Comb2, WSLGs[[x]])))
}

WSOuts <- list()
for(i in 1:12){
  WSOuts[[i]] <- as.data.frame(subset(WSLGs[[i]], Cook>(4/length(WSLGs[[i]]$Cook))))
}
#get all problematic markers
WSPM <- do.call("rbind", WSOuts)
#68

#combine into single file
probs <- data.frame("SNPID"=c(RADPM[,"SNPID"], LM1PM[,"SNPID"], LM2PM[,"SNPID"], CPM[,"SNPID"], WSPM[,"SNPID"]))
#1129

#How many unique SNPs are causing issues?
length(unique(probs$SNPID))
#907
uni4n <- data.frame("SNPID"=unique(probs$SNPID))

#How many are causing repeated problems
probs$ct <- c(rep(1, nrow(probs)))
agg <- aggregate(probs$ct, by=list(probs$SNPID), length)
nrow(agg[which(agg$x>2),])
#29
nrow(agg[which(agg$x>1),])
#191

##Decided to use all markers causing misalignments
write.csv(uni4n, paste(datafiles, file="All_Prob_Cook4n.csv", sep=""))

##STEP 4: Combine problems from 3 sources (gaps, miss-grouped LGs, misaligned positions)
#Find unique list
#extract from list of mapped markers to get list of markers to use in final mapping 

#load data
uni4n <- read.csv(paste(datafiles, "All_Prob_Cook4n.csv", sep=""))
gap_markers <- read.csv(paste(datafiles, "RADSNP_SCLM2_GapMarkers.csv", sep=""))
smallmiss <- read.csv(paste(datafiles, "MissLG_multi.csv", sep=""))
allprobcook <- data.frame("SNPID"=c(uni4n$SNPID, gap_markers$SNPID, smallmiss$Group.1))
#984
uniquecook <- data.frame("SNPID"=c(unique(allprobcook$SNPID)))
#934
write.csv(uniquecook, paste(datafiles, "Problem_withRAD_Cook4n.csv", sep=""))

#Now we want a file of which SNPs to include to subset a vcf file fo rmapping
#bring in all mapped SNPs
comb2 <- read.csv(paste(datafiles, "RADSNP_SCLM2_Map_SNPID.csv", sep=""), header=T)
#subtract those in the problem dataset
toinclude <- data.frame("SNPID"=comb2[-which(comb2$SNPID %in% uniquecook$SNPID),]$SNPID)
#21571

write.table(toinclude, file=paste(datafiles, "RADSNP_SCLM2_RemovedProb_Cook4n.txt", sep=""), col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

#subset using vcftools (see bash file)
