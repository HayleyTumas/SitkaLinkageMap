#Compare Sitka spruce map to each species map
#to get stats on the % of markers that group the same
#and the correlation of marker positions on LGs that group the same

datafiles <- #Path name to DataFiles
output <- #Path to output

#load output from RADSeq to GCAT BLAST and mapped GCAT to species BLASTs
gcat <- read.csv(paste(datafiles, "GCAT_positions_RADChip.csv", sep=""))
limber <- read.csv(paste(datafiles, "Limber_Sitka_p90e100_compare.csv", sep=""))
norway <- read.csv(paste(datafiles, "NorwaySpruce_Sitka_map_compare.csv", sep=""))
#The White spruce map is from Pavy et al. 2017
#No BLAST was needed here because all markers on this map are also based on a GCAT sequence
#allowing for a 1:1 comparison following the RADSeq to GCAT BLAST
ws <- read.csv(paste(datafiles, "White_Spruce.csv", sep=""))

#Rename columns to identify species when merging files
#NOTE: PF=Pinus flexilis (limber pine), PA=Picea abies (Norway spruce), and PS=Picea sitchensis (Sitka spruce)
colnames(limber)[9:10] <- c("LG.PF", "Position.PF")
colnames(norway)[9] <- c("LG.PA")
colnames(norway)[11]<- c("Position.PA")
colnames(gcat)[3:4] <- c("LG.PS", "Position.PS")

##Sitka Spruce to Limber Pine
#Linkage group in limber pine map is a number, so make a factor
limber$LG.PF <- as.factor(limber$LG.PF)
#Determine which linkage group labels match up by aggregating by Sitka and Limber pine linkage groups
#i.e. if 200 of the same markers are on LG1 in Sitka and LG4 in Limber pine, these are essentially the same linkage group
limber_lg <- aggregate(limber$Position.PF, by = limber[c('LG.PS', 'LG.PF')], length)
#subset the full comparison file to include only those markers that group into the same linkage group using the above dataframe 
complimber <- limber[which((limber$LG.PS=="LG1"&limber$LG.PF=="2")|(limber$LG.PS=="LG4"&limber$LG.PF=="1")|(limber$LG.PS=="LG2"&limber$LG.PF=="4")|(limber$LG.PS=="LG6"&limber$LG.PF=="10")|(limber$LG.PS=="LG7"&limber$LG.PF=="3")|(limber$LG.PS=="LG8"&limber$LG.PF=="8")|(limber$LG.PS=="LG3"&limber$LG.PF=="5")|(limber$LG.PS=="LG10"&limber$LG.PF=="9")|(limber$LG.PS=="LG9"&limber$LG.PF=="12")|(limber$LG.PS=="LG5"&limber$LG.PF=="7")|(limber$LG.PS=="LG12"&limber$LG.PF=="6")|(limber$LG.PS=="LG11"&limber$LG.PF=="11")),]
#Check to make sure the subset worked - there should now be only 12 rows in this data frame corresponding to 12 linkage groups
#And each linkge group should have >10 markers (usually >100 markers depending on number of matches in comparison)
limberchk <- aggregate(complimber$Position.PF, by = complimber[c('LG.PS', 'LG.PF')], length)
write.csv(complimber, file=paste(output, "Limber_Match.csv", sep=""))
#Check the % of markers that group on the same chromosome
nrow(complimber)/nrow(limber)
#0.8480845

#Check correlation in marker position for markers that group to the same linkage group
LimberLGs <- split(complimber , f = complimber$LG.PS )
Limbercorr <- list()

for(x in 1:12){
  nam <- paste("LG", x, sep="")
  Limbercorr[[nam]]  <- data.frame("Method"=cor.test(LimberLGs[[x]]$Position.PS, LimberLGs[[x]]$Position.PF, method="kendall", ci=FALSE)$method, "Estimate"=cor.test(LimberLGs[[x]]$Position.PS, LimberLGs[[x]]$Position.PF, method="kendall", ci=FALSE)$estimate, "Pvalue"=as.numeric(cor.test(LimberLGs[[x]]$Position.PS, LimberLGs[[x]]$Position.PF, method="kendall", ci=FALSE)$p.value))
}

Limbercorrs <- do.call("rbind", Limbercorr)
write.csv(Limbercorrs, file=paste(output, "Limber_Correlations.csv", sep=""))
mean(abs(Limbercorrs$Estimate))
#0.9308746
max(abs(Limbercorrs$Estimate))
# 0.9746754
min(abs(Limbercorrs$Estimate))
# 0.8836974

##Norwary Spruce - see Limber Pine code for comments on each step
norway$LG.PA <- as.factor(norway$LG.PA)
norway_lg <- aggregate(norway$Position.PA, by = norway[c('LG.PS', 'LG.PA')], length)
compnorway <- norway[which((norway$LG.PS=="LG1"&norway$LG.PA=="1")|(norway$LG.PS=="LG7"&norway$LG.PA=="3")|(norway$LG.PS=="LG5"&norway$LG.PA=="4")|(norway$LG.PS=="LG4"&norway$LG.PA=="5")|(norway$LG.PS=="LG2"&norway$LG.PA=="2")|(norway$LG.PS=="LG10"&norway$LG.PA=="6")|(norway$LG.PS=="LG6"&norway$LG.PA=="12")|(norway$LG.PS=="LG8"&norway$LG.PA=="8")|(norway$LG.PS=="LG9"&norway$LG.PA=="9")|(norway$LG.PS=="LG3"&norway$LG.PA=="7")|(norway$LG.PS=="LG12"&norway$LG.PA=="10")|(norway$LG.PS=="LG11"&norway$LG.PA=="11")),]
norwaychk <- aggregate(compnorway$Position.PA, by = compnorway[c('LG.PS', 'LG.PA')], length)
write.csv(compnorway, file=paste(output, "Norway_Match.csv", sep=""))
#% of matched markers that are assigned to the same linkage group
nrow(compnorway)/nrow(norway)
#0.8834261

#Some Norway spruce genome sequences are represented multiple times
#The issue is that some SNPs on the same genome sequence were assigned to diffrent linkage groups in the Norwar spruce map
#So will need to remove these duplicate sequences to run statistics
#STATs on matches:
#2857 matches
#1875 unique SNPs
#1755 unique GCAT - NOTE: not all unique GCAT sequences becasue some GCATs have multiple SNPs as well
compnorway_nodups <- compnorway[-which(duplicated(compnorway$SNPID)),] #1875
#do the same for original dataset
#3234 matches
norway_nodups <- norway[-which(duplicated(norway$SNPID)),]
#2071 unique matches
#Check the number of markers that are now assigned to the same linkage group
nrow(compnorway_nodups)/nrow(norway_nodups)
#0.9053597

#check correlation in markers on the same linkage group - using the comparison without duplicates
NorwayLGs <- split(compnorway_nodups , f = compnorway_nodups$LG.PS )
Norwaycorr <- list()

for(x in 1:12){
  nam <- paste("LG", x, sep="")
  Norwaycorr[[nam]]  <- data.frame("Method"=cor.test(NorwayLGs[[x]]$Position.PS, NorwayLGs[[x]]$Position.PA, method="kendall", ci=FALSE)$method, "Estimate"=cor.test(NorwayLGs[[x]]$Position.PS, NorwayLGs[[x]]$Position.PA, method="kendall", ci=FALSE)$estimate, "Pvalue"=as.numeric(cor.test(NorwayLGs[[x]]$Position.PS, NorwayLGs[[x]]$Position.PA, method="kendall", ci=FALSE)$p.value))
}

Norwaycorrs <- do.call("rbind", Norwaycorr)
write.csv(Norwaycorrs, file=paste(output, "Norway_Correlations.csv", sep=""))
mean(abs(Norwaycorrs$Estimate))
# 0.9609283
max(abs(Norwaycorrs$Estimate))
# 0.9823256
min(abs(Norwaycorrs$Estimate))
#0.9110894

##White spruce
#See Limber pine comparison above for comments on code
#merge 2 map files
white <- merge(gcat, ws, by="GCAT") #NOTE: If multiple SNPs on a GCAT, both SNPs included, i.e. duplicated GCATs
write.csv(white, file=paste(output, "White_Sitka_map_compare.csv", sep="")) 
white_lg <- aggregate(white$Position.WS, by = white[c('LG.PS', 'LG.WS')], length)
compwhite <- white[which((white$LG.PS=="LG1"&white$LG.WS=="LG07")|(white$LG.PS=="LG2"&white$LG.WS=="LG08")|(white$LG.PS=="LG4"&white$LG.WS=="LG01")|(white$LG.PS=="LG3"&white$LG.WS=="LG02")|(white$LG.PS=="LG6"&white$LG.WS=="LG05")|(white$LG.PS=="LG7"&white$LG.WS=="LG11")|(white$LG.PS=="LG8"&white$LG.WS=="LG06")|(white$LG.PS=="LG5"&white$LG.WS=="LG03")|(white$LG.PS=="LG9"&white$LG.WS=="LG04")|(white$LG.PS=="LG10"&white$LG.WS=="LG12")|(white$LG.PS=="LG12"&white$LG.WS=="LG09")|(white$LG.PS=="LG11"&white$LG.WS=="LG10")),]
whitechk <- aggregate(compwhite$Position.WS, by = compwhite[c('LG.PS', 'LG.WS')], length)
write.csv(compwhite, file=paste(output, "White_Match.csv", sep=""))
#% markers on the same linkage group
nrow(compwhite)/nrow(white)
#0.9359251
#check correlation in markers grouped on the same linkage group
WhiteLGs <- split(compwhite , f = compwhite$LG.PS )
Whitecorr <- list()

for(x in 1:12){
  nam <- paste("LG", x, sep="")
  Whitecorr[[nam]]  <- data.frame("Method"=cor.test(WhiteLGs[[x]]$Position.PS, WhiteLGs[[x]]$Position.WS, method="kendall", ci=FALSE)$method, "Estimate"=cor.test(WhiteLGs[[x]]$Position.PS, WhiteLGs[[x]]$Position.WS, method="kendall", ci=FALSE)$estimate, "Pvalue"=as.numeric(cor.test(WhiteLGs[[x]]$Position.PS, WhiteLGs[[x]]$Position.WS, method="kendall", ci=FALSE)$p.value))
}

Whitecorrs <- do.call("rbind", Whitecorr)
write.csv(Whitecorrs, file=paste(output, "White_Correlations.csv", sep=""))
mean(abs(Whitecorrs$Estimate))
#  0.9758238
max(abs(Whitecorrs$Estimate))
# 0.9833165
min(abs(Whitecorrs$Estimate))
#0.9590509

#Create map of comparison between Sitka and white spruce
#NOTE: graphs for other species are in output, code for these graphs are in species BLAST folders
library(ggplot2)
library(ggpubr)

white$LG.PS <- factor(white$LG.PS, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
PS_PG <- ggplot(white, aes(x=Position.PS, y=Position.WS))+geom_point(size=1)+
  ylab("Position in White Spruce")+ xlab("Position in Sitka Spruce")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_grid(LG.WS~LG.PS)


#Check for inversions when comparing maps between species - looking at markers on the same linkage group
limbermatches <- read.csv(paste(datafiles, "Limber_Match.csv", sep=""))
norwaymatches <- read.csv(paste(datafiles, "Norway_Match.csv", sep="")) 
whitematches <- read.csv(paste(datafiles, "White_Match.csv", sep=""))

#Limber Pine
#Need to order dataframe based on Sitka LG and position to be able to compare position
limbermatches <- limbermatches[order(limbermatches$LG.PS, limbermatches$Position.PS),]
#Assign SNP number rather than position that can be used to order SNPs based on Sitka position
limbermatches$SNPNUM <- c(1:nrow(limbermatches))
#Break out by species to format dataframe for graphing
limber_pf <- limbermatches[, c("SNPID", "SNPNUM", "LG.PS", "Position.PF")]
limber_ps <- limbermatches[, c("SNPID", "SNPNUM", "LG.PS", "Position.PS")]
colnames(limber_pf)[4] <- "Position"
colnames(limber_ps)[4] <- "Position"
#Add column identifying species
limber_pf$Species <- "P.Flexilis"
limber_ps$Species <- "P.sitchensis"
#Combine species datasets
limbercomp <- rbind(limber_pf, limber_ps)
#Order linkge groups in correct order
limbercomp$LG.PS <- factor(limbercomp$LG.PS, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))

limberin <- ggplot(limbercomp, aes(x=SNPNUM, y=Species, color=Position))+geom_point(size=1)+
  ylab("Species")+ xlab("Position")+
  theme_bw()+scale_color_gradientn(colors=rainbow(5))+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+theme(legend.position = "none")+
  facet_wrap(~LG.PS, scales="free")

ggsave(paste(output,"Sitka_Limber_Inversion_Chk.png", sep=""), limberin, height=5, width=10, units="in")


#Norway Spruce
norwaymatches <- norwaymatches[order(norwaymatches$LG.PS, norwaymatches$Position.PS),]
norwaymatches$SNPNUM <- c(1:nrow(norwaymatches))
norway_pa <- norwaymatches[, c("SNPID", "SNPNUM", "LG.PS", "Position.PA")]
norway_ps <- norwaymatches[, c("SNPID", "SNPNUM", "LG.PS", "Position.PS")]
colnames(norway_pa)[4] <- "Position"
colnames(norway_ps)[4] <- "Position"
norway_pa$Species <- "P.abies"
norway_ps$Species <- "P.sitchensis"
norwaycomp <- rbind(norway_pa, norway_ps)
norwaycomp$LG.PS <- factor(norwaycomp$LG.PS, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
  
norwayin <- ggplot(norwaycomp, aes(x=SNPNUM, y=Species, color=Position))+geom_point(size=1)+
  ylab("Species")+ xlab("Position")+
  theme_bw()+scale_color_gradientn(colors=rainbow(5))+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+theme(legend.position = "none")+
  facet_wrap(~LG.PS, scales="free")

ggsave(paste(output,"Sitka_Norway_Inversion_Chk.png", sep=""), norwayin, height=5, width=10, units="in")

##White Spruce
whitematches <- whitematches[order(whitematches$LG.PS, whitematches$Position.PS),]
whitematches$SNPNUM <- c(1:nrow(whitematches))
white_pg <- whitematches[, c("SNPID.x", "SNPNUM", "LG.PS", "Position.WS")]
white_ps <- whitematches[, c("SNPID.x", "SNPNUM", "LG.PS", "Position.PS")]
colnames(white_pg)[4] <- "Position"
colnames(white_ps)[4] <- "Position"
white_pg$Species <- "P.glauca"
white_ps$Species <- "P.sitchensis"
whitecomp <- rbind(white_pg, white_ps)
whitecomp$LG.PS <- factor(whitecomp$LG.PS, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))

whitein <- ggplot(whitecomp, aes(x=SNPNUM, y=Species, color=Position))+geom_point(size=1)+
  ylab("Species")+ xlab("Position")+
  theme_bw()+scale_color_gradientn(colors=rainbow(5))+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+theme(legend.position = "none")+
  facet_wrap(~LG.PS, scales="free")

ggsave(paste(output, "Sitka_White_Inversion_Chk.png", sep=""), whitein, height=5, width=10, units="in")
