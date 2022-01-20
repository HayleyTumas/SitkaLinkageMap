#Need to compare map from RADSeq and Chip data with problematic markers removed
#to the 4 component maps: maps from SNP Chip and RADSeq separately
#and maps from SNP Chip and RADseq data together for each mapping family
#Will get stats on the % of markers that group the same
#and the correlation of marker positions on LGs that group the same

datafiles <- #Path to DataFiles
output <- #Path to Output

#load data
#Map from RADSeq and Chip data with problematic markers removed (Full map):
comb <- read.csv(paste(datafiles, "RADSNP_SCLM2_Cook4NARC_Map_SNPID.csv", sep=""), header=T)
#load other component maps
#Map using RADSeq and Chip data for only Family 1:
lm1 <- read.csv(paste(datafiles, "RADSNP_LM1_Map_SNPID.csv", sep=""), header=T)
#Map using RADSeq and Chip data for only Family 2:
lm2 <- read.csv(paste(datafiles, "RADSNP_LM2_Map_SNPID.csv", sep=""), header=T)
#Map using both families, but only SNP Chip data:
chip <- read.csv(paste(datafiles, "CombQCNR_LM2SC_Map_SNPID.csv", sep=""))
#Map using both families, but only RADSeq data:
RAD <- read.csv(paste(datafiles, "RADQC_2FAMSSSCLM2_Map_SNPID.csv", sep=""))
#NOTE: Code for creating maps can be found in respective folders

#Rename columns to identify maps when merged
colnames(lm1)[2:5] <- paste(colnames(lm1)[2:5], "LM1", sep=".")
colnames(lm2)[2:5] <- paste(colnames(lm2)[2:5], "LM2", sep=".")
colnames(comb)[2:5] <- paste(colnames(comb)[2:5], "Comb", sep=".")
colnames(chip)[2:5] <- paste(colnames(chip)[2:5], "Chip", sep=".")
colnames(RAD)[2:5] <- paste(colnames(RAD)[2:5], "RAD", sep=".")

##Full map:Family 1 Map
#combine files
lm1comb <- merge(comb, lm1, by="SNPID")
#aggregate data by linkage group in both maps
#this helps identify which linkage groups match in each map
#i.e. if 300 of the same markers are on LG1 in the full map and LG5 in family 1 map
#then LG1 and LG5 are the same chromosome 
#Some markers are assigned to different linkage groups
#this would appear as 8 markers being assigned to LG1 in the full map and LG8 in family 1, for example
lm1comb_lg <- aggregate(lm1comb$Position.LM1, by = lm1comb[c('LG.Comb', 'LG.LM1')], length)
#subset data based on markers that group on the same liknage group - using above dataframe
complm1comb <- lm1comb[which((lm1comb$LG.Comb=="LG1"&lm1comb$LG.LM1=="LG1")|(lm1comb$LG.Comb=="LG10"&lm1comb$LG.LM1=="LG9")|(lm1comb$LG.Comb=="LG11"&lm1comb$LG.LM1=="LG12")|(lm1comb$LG.Comb=="LG12"&lm1comb$LG.LM1=="LG11")|(lm1comb$LG.Comb=="LG2"&lm1comb$LG.LM1=="LG8")|(lm1comb$LG.Comb=="LG3"&lm1comb$LG.LM1=="LG3")|(lm1comb$LG.Comb=="LG4"&lm1comb$LG.LM1=="LG2")|(lm1comb$LG.Comb=="LG5"&lm1comb$LG.LM1=="LG7")|(lm1comb$LG.Comb=="LG6"&lm1comb$LG.LM1=="LG4")|(lm1comb$LG.Comb=="LG7"&lm1comb$LG.LM1=="LG6")|(lm1comb$LG.Comb=="LG8"&lm1comb$LG.LM1=="LG10")|(lm1comb$LG.Comb=="LG9"&lm1comb$LG.LM1=="LG5")),]
#check the subset work - there should be only 12 rows in this new dataframe
#representing the 12 chromosomes (regardless of linkage group names in each map)
#And each chromosome should have >100 markers
lm1chk <- aggregate(complm1comb$Position.LM1, by = complm1comb[c('LG.Comb', 'LG.LM1')], length)
write.csv(complm1comb, file="LM1_Match.csv")
#% of markers assigned to the same linkage group
nrow(complm1comb)/nrow(lm1comb)
#0.9986011
#check correlation in marker position for markers that group to the same linkage group
LM1LGs <- split(complm1comb , f = complm1comb$LG.Comb )
LM1corr <- list()

for(x in 1:12){
  nam <- paste("LG", x, sep="")
  LM1corr[[nam]]  <- data.frame("Method"=cor.test(LM1LGs[[x]]$Position.Comb, LM1LGs[[x]]$Position.LM1, method="kendall", ci=FALSE)$method, "Estimate"=cor.test(LM1LGs[[x]]$Position.Comb, LM1LGs[[x]]$Position.LM1, method="kendall", ci=FALSE)$estimate, "Pvalue"=as.numeric(cor.test(LM1LGs[[x]]$Position.Comb, LM1LGs[[x]]$Position.LM1, method="kendall", ci=FALSE)$p.value))
}

LM1corrs <- do.call("rbind", LM1corr)
write.csv(LM1corrs, file="LM1_Correlations.csv")
mean(abs(LM1corrs$Estimate))
#0.985638
min(abs(LM1corrs$Estimate))
#0.9758894
max(abs(LM1corrs$Estimate))
# 0.991182

#Full Map:Family 2 map
#See Family 1 comparison for code comments
lm2comb <- merge(comb, lm2, by="SNPID")
lm2comb_lg <- aggregate(lm2comb$Position.LM2, by = lm2comb[c('LG.Comb', 'LG.LM2')], length)
complm2comb <- lm2comb[-which(lm2comb$LG.Comb=="LG6"&lm2comb$LG.LM2=="LG2"),]
lm2chk <- aggregate(complm2comb$Position.LM2, by = complm2comb[c('LG.Comb', 'LG.LM2')], length)
write.csv(complm2comb, file="LM2_Match.csv")
nrow(complm2comb)/nrow(lm2comb)
#0.9999298
#check correlation in position of markers that group on the same linkage group
LM2LGs <- split(complm2comb , f = complm2comb$LG.Comb )
LM2corr <- list()

for(x in 1:12){
  nam <- paste("LG", x, sep="")
  LM2corr[[nam]]  <- data.frame("Method"=cor.test(LM2LGs[[x]]$Position.Comb, LM2LGs[[x]]$Position.LM2, method="kendall", ci=FALSE)$method, "Estimate"=cor.test(LM2LGs[[x]]$Position.Comb, LM2LGs[[x]]$Position.LM2, method="kendall", ci=FALSE)$estimate, "Pvalue"=as.numeric(cor.test(LM2LGs[[x]]$Position.Comb, LM2LGs[[x]]$Position.LM2, method="kendall", ci=FALSE)$p.value))
}

LM2corrs <- do.call("rbind", LM2corr)

write.csv(LM2corrs, file="LM2_Correlations.csv")
mean(abs(LM2corrs$Estimate))
#0.9770212
min(abs(LM2corrs$Estimate))
#0.9572067
max(abs(LM2corrs$Estimate))
#0.990487

#Full map:Map from SNP Chip data
chipcomb <- merge(comb, chip, by="SNPID")
chipcomb_lg <- aggregate(chipcomb$Position.Chip, by = chipcomb[c('LG.Comb', 'LG.Chip')], length)
compchipcomb <- chipcomb[-which(chipcomb$LG.Comb=="LG6"&chipcomb$LG.Chip=="LG8"),]
chipchk <- aggregate(compchipcomb$Position.Chip, by = compchipcomb[c('LG.Comb', 'LG.Chip')], length)
write.csv(compchipcomb, file="Chip_Match.csv")
#% lg match
nrow(compchipcomb)/nrow(chipcomb)
#0.999781
#check correlation in LGs that match
ChipLGs <- split(compchipcomb , f = compchipcomb$LG.Comb )
Chipcorr <- list()

for(x in 1:12){
  nam <- paste("LG", x, sep="")
  Chipcorr[[nam]]  <- data.frame("Method"=cor.test(ChipLGs[[x]]$Position.Comb, ChipLGs[[x]]$Position.Chip, method="kendall", ci=FALSE)$method, "Estimate"=cor.test(ChipLGs[[x]]$Position.Comb, ChipLGs[[x]]$Position.Chip, method="kendall", ci=FALSE)$estimate, "Pvalue"=as.numeric(cor.test(ChipLGs[[x]]$Position.Comb, ChipLGs[[x]]$Position.Chip, method="kendall", ci=FALSE)$p.value))
}

Chipcorrs <- do.call("rbind", Chipcorr)
write.csv(Chipcorrs, file="Chip_Correlations.csv")
mean(abs(Chipcorrs$Estimate))
#0.9500679
min(abs(Chipcorrs$Estimate))
#0.8704772 - LG6
max(abs(Chipcorrs$Estimate))
# 0.9859884

#Full Map:Map using only RADSeq data
radcomb <- merge(comb, RAD, by="SNPID")
radcomb_lg <- aggregate(radcomb$Position.RAD, by = radcomb[c('LG.Comb', 'LG.RAD')], length)
compradcomb <- radcomb[-which(radcomb$LG.Comb=="LG6"&radcomb$LG.RAD=="LG7"),]
radchk <- aggregate(compradcomb$Position.RAD, by = compradcomb[c('LG.Comb', 'LG.RAD')], length)
write.csv(compradcomb, file="RAD_Match.csv")
nrow(compradcomb)/nrow(radcomb)
# 0.9995817
#check correlation in position of markers grouped on same linkage group
RADLGs <- split(compradcomb , f = compradcomb$LG.Comb )
RADcorr <- list()

for(x in 1:12){
  nam <- paste("LG", x, sep="")
  RADcorr[[nam]]  <- data.frame("Method"=cor.test(RADLGs[[x]]$Position.Comb, RADLGs[[x]]$Position.RAD, method="kendall", ci=FALSE)$method, "Estimate"=cor.test(RADLGs[[x]]$Position.Comb, RADLGs[[x]]$Position.RAD, method="kendall", ci=FALSE)$estimate, "Pvalue"=as.numeric(cor.test(RADLGs[[x]]$Position.Comb, RADLGs[[x]]$Position.RAD, method="kendall", ci=FALSE)$p.value))
}

RADcorrs <- do.call("rbind", RADcorr)
write.csv(RADcorrs, file="RAD_Correlations.csv")
mean(abs(RADcorrs$Estimate))
# 0.9803597
min(abs(RADcorrs$Estimate))
#0.9687906
max(abs(RADcorrs$Estimate))
# 0.9847505


####Figures
library(ggplot2)
library(ggpubr)

lg.levels <- c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12")
LGs <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")

#Full:Family 1
lm1comb$LG.Comb <- factor(lm1comb$LG.Comb, levels=lg.levels)
lm1comb$LG.LM1 <- factor(lm1comb$LG.LM1, levels=lg.levels)
LM1 <- ggplot(lm1comb, aes(x=Position.Comb, y=Position.LM1))+geom_point(size=1)+
  ylab("Position in Mapping Family 1")+ xlab("Position in Composite Map")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_grid(LG.LM1~LG.Comb)

ggsave(paste(outptut, "Full_Family1.png", sep=""), LM1, height=6, width=6, units="in")

#Full:Family 2
lm2comb$LG.Comb <- factor(lm2comb$LG.Comb, levels=lg.levels)
lm2comb$LG.LM2 <- factor(lm2comb$LG.LM2, levels=lg.levels)
LM2 <- ggplot(lm2comb, aes(x=Position.Comb, y=Position.LM2))+geom_point(size=1)+
  ylab("Position in Mapping Family 2")+ xlab("Position in Composite Map")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_grid(LG.LM2~LG.Comb)

ggsave(paste(outptut,"Full_Family2.png", sep=""), LM2, height=6, width=6, units="in")

#Full:Chip
chipcomb$LG.Comb <- factor(chipcomb$LG.Comb, levels=lg.levels)
chipcomb$LG.Chip <- factor(chipcomb$LG.Chip, levels=lg.levels)
ChipGraph <- ggplot(chipcomb, aes(x=Position.Comb, y=Position.Chip))+geom_point(size=1)+
  ylab("Position in SNP Chip Map")+ xlab("Position in Composite Map")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_grid(LG.Chip~LG.Comb)

ggsave(paste(outptut,"Full_Chip.png", sep=""), ChipGraph, height=6, width=6, units="in")

#Full:RADSeq
radcomb$LG.Comb <- factor(radcomb$LG.Comb, levels=lg.levels)
radcomb$LG.RAD <- factor(radcomb$LG.RAD, levels=lg.levels)
RADGraph <- ggplot(radcomb, aes(x=Position.Comb, y=Position.RAD))+geom_point(size=1)+
  ylab("Position in RADSeq Map")+ xlab("Position in Composite Map")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_grid(LG.RAD~LG.Comb)

ggsave(paste(outptut,"Full_RAD.png", sep=""), RADGraph, height=6, width=6, units="in")

