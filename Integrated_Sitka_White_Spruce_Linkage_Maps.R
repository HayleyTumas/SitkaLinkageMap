#One goal of the project was to create a map that integrated markers from the white spruce and Sitka spruce maps.
#Because both maps had markers based on the white spruce gene catalog (GCAT),
#there were a set of markers found on both maps that could be used as reference points
#for integrating the maps into a single consensus map. 
#We used LPMerge, a package that is usually used to create consensus maps within a species
#from multiple population or family component maps.
#In order to integrate the maps using LPMerge, we needed to:
#1)  Select one SNP per gene catalog locus in the Sitka spruce map, preferentially selecting 
#    SNP Chip markers that were linked to gene catalog sequences through exome capture. Then 
#    rename all Sitka spruce markers with a gene catalog match as the gene catalog sequence to be 
#    comparable ot the white spruce map.
#2)  Remove any markers that disrupted synteny between the 2 species maps. 
#    This includes markers that did not group on the same linkage group or misaligned.
#3)  Format the data for LPMerge and order linkage groups in the same order across the 2 species.
#3B) Create input for only markers that were mapped in both species to act as a reference.
#4)  Run LPMerge across multiple intervals on each linkage group for the subset of markers mapped in both species.
#5)  Select the best interval for each linkage group based on RMSE and linkage group length.
#6)  Repeat step 4 and 5 using the integrated map of markers mapped in both species and each species maps as components.
#7)  Manually emove gaps created by the longer linkage groups in Sitka spruce.
#8A) Calculate final integrated map statistics. 
#8B) Test map synteny with the component species maps

#set working directories
datadir <- "paste path to the Data folder" 
outputdir <- "path to output folder"

#Load Data - files comparing maps were created using code from Species_Comparison folder
#File names here are the same file names in the respective code files in the Species_Comparison folder
#file comparing linkage group and position of white and Sitka maps (Compare_Sitka_toSpecies.R)
intmarkers <- read.csv(paste(datadir, "White_Sitka_map_compare.csv", sep="/"))
#file with markers found on the same linkage group between the two maps (Compare_Sitka_toSpecies.R)
sameLG <- read.csv(paste(datadir, "White_Match.csv", sep="/"))
#finalized Sitka spruce RAD-Chip map
sitka <- read.csv(paste(datadir, "RADSNP_SCLM2_Cook4NARC_Map_SNPID.csv", sep="/"))
#white spruce map from Pavy et al. 2017
#this includes the GCAT positions of all SNPs in the SNP ID
white <- read.csv(paste(datadir, "White_Spruce.csv", sep="/"))
#GCAT positions of all mapped Sitka SNPs (both RAD and Chip)
#Note: SNP Chip SNPs in this file were linked to GCAT sequenuces during SNP discovery
#This file that includes RAD Seq SNPs was created with code in the Species_Comparisons folder
#RADSeq_GCAT code
sitka_gcat <- read.csv(paste(datadir, "GCAT_positions_RADChip.csv", sep="/"))

##STEP 1: Get Sitka map with only one SNP per GCAT (so no duplicated GCATs)
#merge Sitka map with Sitka GCAT key
sitka_map_gcat <- merge(sitka, sitka_gcat, by="SNPID", all.x=TRUE) #21570 (all unique)
#now select 1 SNP per GCAT, and preferentially remove RADSeq SNPs from duplicated GCATs
#because SNP Chip SNPs have stronger link to GCAT, whereas RADSeq SNPs are linked based on BLAST
#creat column that designates RADSeq versus Chip SNPs
#All RADSeq SNPs have "SNP" in ID, no Chip SNPs have this in ID
sitka_map_gcat$Type[(grepl("SNP", sitka_map_gcat$SNPID, fixed=TRUE))] <- "RAD"
sitka_map_gcat$Type[is.na(sitka_map_gcat$Type)] <- "CHIP"
#Order dataframe based on this column to remove RADSeq SNPs first 
#(removes second entry and RAD is alphabetically second)
sitka_type <- sitka_map_gcat[order(sitka_map_gcat$Type),]
#Get dataframe of any SNP without a GCAT sequence (all RADSeq SNPs)
#Because GCAT sequences would be NA, removing duplicates would remove all of these SNPs
sitka_noGCAT <- sitka_type[is.na(sitka_type$GCAT),] #15812
#Get only SNPs with a GCAT sequence
sitka_wGCAT <- sitka_type[-which(sitka_type$SNPID %in% sitka_noGCAT$SNPID),] 
#5758, 5414 unique GCAT
#Remove duplicated GCAT sequences
sitka_single <- sitka_wGCAT[-which(duplicated(sitka_wGCAT$GCAT)),]
#344 removed, 5414 remain
#Recombine SNPs without a GCAT and the new file with each GCAT represented by a single SNP
sitka_singleGCAT <- rbind(sitka_single, sitka_noGCAT) 
#21,226 SNPs remain
#Remove any markers that did not group to the same linkage group across the two maps
#Need to do it by SNP ID because for a few GCAT with multiple SNPs, those SNPs were assigned to different linkage groups
#So a GCAT may be in sameLG even if all SNPs associated with the GCAT did not map to the same linkage group
badmatch <- intmarkers[-which(intmarkers$SNPID.x %in% sameLG$SNPID.x),]#178 total, 176 unique GCAT
sitka_goodsingle <- sitka_singleGCAT[-which(sitka_singleGCAT$SNPID %in% badmatch$SNPID.x),] 
#21,061 SNPs remain
#Create a column that has GCAT, if GCAT, or SNPID if no GCAT for integrating with white spruce
#Since the GCAT ID will what be used as scaffolds or reference SNPs in LPMerge
#Will need to have teh same ID in both files
sitka_goodsingle$MarkerID <- ifelse(is.na(sitka_goodsingle$GCAT), sitka_goodsingle$SNPID, sitka_goodsingle$GCAT) #21,009 unique marker ids

#STEP 2: Eliminate markers that did not align between the two maps
#Using only markers that group on the same linkage group from the Species_Comparison folder 
#and removing any markers that have a different position on the same linkage group
#based on Cook's distance (as when identifying problematic markers in the RAD-Chip map)
white_rmbad <- white[-which(white$GCAT %in% badmatch$GCAT),]#8617 (see badmatch in STEP 1)
sitka_white_input <- merge(sitka_goodsingle, white_rmbad, by="GCAT")#2405
#Check to find those on the same linkage group
input_agg <- aggregate(sitka_white_input$GCAT, by=sitka_white_input[,c("LG.x", "LG.WS")], length)
#all group the same due to filtering from comparison
#Now need to find any markers that substantially misalign based on Cook's distance
LGs <- split(sitka_white_input, f = sitka_white_input$LG.x)
for(x in 1:12){
  LGs[[x]]$Cook  <- (cooks.distance(lm(Position.WS~Position.x, LGs[[x]])))
}
Outs <- list()
for(i in 1:12){
  Outs[[i]] <- as.data.frame(subset(LGs[[i]], Cook>(4/length(LGs[[i]]$Cook))))
}
#get all problematic markers
Misalign <- do.call("rbind", Outs)#78
#Remove these markers to create input files
sitka_int <- sitka_goodsingle[-which(sitka_goodsingle$SNPID %in% Misalign$SNPID.x),]#20983 (removed 78)
white_int <- white_rmbad[-which(white_rmbad$GCAT %in% Misalign$GCAT),]#8539 (removed 78)
#check overlap
inputoverlap <- merge(sitka_int, white_int, by="GCAT") #2327
#check that it is only unique GCATs
length(unique(inputoverlap$GCAT))
#2327
#check group the same
groupsame <- aggregate(inputoverlap$GCAT, by=list(inputoverlap$LG.x, inputoverlap$LG.WS), length)
#all group the same

##STEP 3: Format map files for LPMerge
#Need to have a file with only LG, position, and Marker ID
sitka_lp <- data.frame("LG"=sitka_int$LG.x, "Marker"=sitka_int$MarkerID, "Position"=sitka_int$Position.x)
#NOTE: GCAT ID rather than SNPID will be used as Maker ID for White Spruce for comparison to Sitka
white_lp <- data.frame("LG"=white_int$LG.WS, "Marker"=white_int$GCAT, "Position"=white_int$Position.WS)
#Save these as input files
write.csv(sitka_lp, file=paste(datadir, "Sitka_LPMerge_Input.csv", sep="/"), quote = F, row.names = F)
write.csv(white_lp, file=paste(datadir, "White_LPMerge_Input.csv", sep="/"), quote= F, row.names = F)

##STEP 3B: Create input of only markers that are found on both maps
overlap <- merge(sitka_lp, white_lp, by="Marker")
sitka_lp_overlap <- sitka_lp[which(sitka_lp$Marker %in% overlap$Marker),]
white_lp_overlap <- white_lp[which(white_lp$Marker %in% overlap$Marker),]
#Save these as input files
write.csv(sitka_lp_overlap, file=paste(datadir, "Sitka_LPMerge_Input_Overlap.csv", sep="/"), quote = F, row.names = F)
write.csv(white_lp_overlap, file=paste(datadir, "White_LPMerge_Input_Overlap.csv", sep="/"), quote= F, row.names = F)

#Need to get a separate dataframe for each linkage group
SitkaLGs <- split( sitka_lp_overlap , f = sitka_lp_overlap$LG)
WhiteLGs <- split( white_lp_overlap , f = white_lp_overlap$LG)

#Then need to remove linkage group column
sitka_input <- lapply(SitkaLGs, function(x) { x["LG"] <- NULL; x })
white_input <- lapply(WhiteLGs, function(x) { x["LG"] <- NULL; x })

#Need to put linkage groups in the same order for each species
#i.e. if Sitka LG1 is the same as White LG4, the both need to be first in the linkage group list
#this will simplify merging each linkage group
#Aggregate table as key:
#     LG.PS LG.WS   x
# 7     LG4  LG01 240
# 17    LG3  LG02 239
# 26    LG5  LG03 203
# 35    LG9  LG04 198
# 43    LG6  LG05 233
# 56    LG8  LG06 205
# 57    LG1  LG07 324
# 70    LG2  LG08 241
# 78   LG12  LG09 168
# 85   LG11  LG10 141
# 101   LG7  LG11 211
# 104  LG10  LG12 197

#Get white spruce linkage groups in matching order
white_input_order <- list(white_input[["LG07"]], white_input[["LG08"]], white_input[["LG02"]], white_input[["LG01"]], white_input[["LG03"]], white_input[["LG05"]], white_input[["LG11"]], white_input[["LG06"]], white_input[["LG04"]], white_input[["LG12"]], white_input[["LG10"]], white_input[["LG09"]])

#AND because R puts numbered linkage groups in numerical order (reading 1:9 incorrectly because there is no zero)
#Fix sitka input order
sitka_input_order <- list(sitka_input[["LG1"]], sitka_input[["LG2"]], sitka_input[["LG3"]], sitka_input[["LG4"]], sitka_input[["LG5"]], sitka_input[["LG6"]], sitka_input[["LG7"]], sitka_input[["LG8"]], sitka_input[["LG9"]], sitka_input[["LG10"]], sitka_input[["LG11"]], sitka_input[["LG12"]])

##Check that all data has been formatted properly
#Check that the Sitka input has only unique markers
#(proving that renaming SNPs and removing duplicate SNPs on GCATs worked correctly)
rows <- list()
unique <- list()
for(x in 1:12){
  rows[[x]] <- data.frame("LG"=paste(x), "Num"=nrow(sitka_input[[x]]))
  unique[[x]] <- data.frame("LG"=paste(x), "Num"=length(unique(sitka_input[[x]][,1])))
}
rowtot <- do.call(rbind, rows)
unitot <- do.call(rbind, unique)
check_marks <- cbind(rowtot, unitot)
colnames(check_marks) <- c("LG", "Num.Row", "LG", "Num.Uni")
length(check_marks$Num.Row==check_marks$Num.Uni)
#Should be 12
sum(unitot$Num)
#should be 2327

#Check that the ordered LG lists are in the correct order
#The total number of overlapping markers when this is checked by linkage group
#should be equal to the total number of overlapping markers in the two input files (2327)
#check
chk <- list()
for(g in 1:12){
  chk[[g]] <- data.frame("LG"=paste(g), "Overlap"=nrow(sitka_input_order[[g]][which(sitka_input_order[[g]][,1] %in% white_input_order[[g]][,1]),]))
}
chktot <- do.call(rbind, chk)
sum(chktot$Overlap)
#Should be 2327

##STEP 4: Run overlapping markers in LPMerge to integrate maps by linkage group
install.packages("LPmerge")
library(LPmerge)

#Chromosome integration is run multiple times across an interval range
#The best integration is based on the lowest RMSE - however low RMSE is not always acheived
#In this case, both RMSE and chromosome length should be used to select the best interval
#With chromosome length being the most similar to the 2 input chromosomes
#Unfortunately the RMSE output cannot be saved as output by LPMerge
#So each linkage group must be run separately and the RMSE output must be saved manually
#When integrating, each map can also be weighted when determining marker order
#Here we weighted the two species maps by samples size used to make the maps (528 for Sitka, 1976 for White)
#For loop - writing all output to a text file
sink(file=paste(outputdir, "LPMerge_Output_Overlap_weighted.txt", sep="/"))
outlistw <- list()
for(x in 1:12){
  print(paste("------------------------","CHR", x, "----------------------------", sep=""))
  nam <- paste("CHR", x, sep="")
  outlistw[[nam]] <- LPmerge(list(sitka_input_order[[x]], white_input_order[[x]]), max.interval=1:10, weights=c(528, 1976))
}
sink()

#visualize output using graphs
outw <- outlistw
chrsw <- list()
for(n in 1:12){
  for(x in 1:10){
    outw[[n]][[x]]$Interval <- rep(x, nrow(outw[[n]][[x]]))
  }
  chrsw[[n]] <- do.call("rbind", outw[[n]])
  chrsw[[n]]$Interval <- as.factor(chrsw[[n]]$Interval)
  print(ggplot(chrsw[[n]], aes(x=consensus, y=Interval))+geom_point()+labs(title=paste("CHR", n))+geom_vline(xintercept=c(maxlength[n,4], max(chrsw[[n]][,3]), max(chrsw[[n]][,4])), color=c("black", "red", "blue")))
  ggsave(paste(paste(outputdir, "/overlap_plots/", sep=""), "CHR", n, "_weighted.png", sep=""), height = 10, width=20, units="cm")
  chrsw[[n]]$LG <- paste("LG", n, sep="")
}

#save full output with all intervals
chrsw_all <- do.call("rbind", chrsw)
write.csv(chrsw_all, paste(outputdir, "overlap_markers_output_weighted.csv", sep="/"), row.names = F, quote=F)

##STEP 5: Select the best interval for each linkage group based on RMSE and linkage group length.
#reload full map output
ovout <- read.csv(paste(outputdir, "overlap_markers_output_weighted.csv", sep="/"))
#select best intervals (based on length, RMSE, or both) based on text file output
##CHR1=10 (length)
##CHR2=4 (RMSE)
##CHR3=10 (length & RMSE)
##CHR4=9 (length & RMSE)
##CHR5=9 (length & RMSE)
##CHR6=7 (length & RMSE)
##CHR7=10 (length & RMSE)
##CHR8=9 (RMSE)
##CHR9=10 (length & RMSE)
##CHR10=2 (RMSE)
##CHR11=2 (RMSE)
##CHR12=10 (length & RMSE)

#create dataframe with only the results for the best interval for each chromosome
ovout_best <- ovout[which(ovout$LG=="LG1"&ovout$Interval=="10"|ovout$LG=="LG2"&ovout$Interval=="4"|
  ovout$LG=="LG3"&ovout$Interval=="10"|ovout$LG=="LG4"&ovout$Interval=="9"|
  ovout$LG=="LG5"&ovout$Interval=="9"|ovout$LG=="LG6"&ovout$Interval=="7"|
    ovout$LG=="LG7"&ovout$Interval=="10"|ovout$LG=="LG8"&ovout$Interval=="9"|
    ovout$LG=="LG9"&ovout$Interval=="10"|ovout$LG=="LG10"&ovout$Interval=="2"|
    ovout$LG=="LG11"&ovout$Interval=="2"|ovout$LG=="LG12"&ovout$Interval=="10"),]

#save this file as the finalized overlap integrated map
write.csv(ovout_best, file=paste(outputdir, "overlap_markers_weighted_bestchrs.csv", sep="/"))

#visually compare overlap map to each species component map using graphs
#to ensure integration is reasonable
library(tidyr)
colnames(ovout_best)[2:4] <- c("consensus", "sitka", "white")
ovout_long <- gather(ovout_best, "Species", "Position", consensus:white)
ovout_long$SpLG <- paste(ovout_long$LG, ovout_long$Species, sep="_")
compare <- ggplot(ovout_long, aes(x=Position, y=SpLG, color=Species))+geom_point()

ovout_best$LG <- factor(ovout_best$LG, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
C_PS <- ggplot(ovout_best, aes(x=consensus, y=sitka))+geom_point(size=1)+
  ylab("Position in Sitka Spruce")+ xlab("Position in Consensus Map")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_wrap(~LG)
C_PG <- ggplot(ovout_best, aes(x=consensus, y=white))+geom_point(size=1)+
  ylab("Position in White Spruce")+ xlab("Position in Consensus Map")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_wrap(~LG)

ggsave(paste(paste(outputdir, "/overlap_plots/", sep=""), "compare_weighted.png", sep=""), compare, height = 20, width=25, units="cm")
ggsave(paste(paste(outputdir, "/overlap_plots/", sep=""), "Sitka_weighted.png", sep=""), C_PS, height = 15, width=20, units="cm")
ggsave(paste(paste(outputdir, "/overlap_plots/", sep=""), "White_weighted.png", sep=""), C_PG, height = 15, width=20, units="cm")

#Format overlap integrated map for input into LPMerge 
combo_lp <- data.frame("LG"=ovout_best$LG, "Marker"=ovout_best$marker, "Position"=ovout_best$consensus)
write.csv(combo_lp, paste(datadir, "Consensus_Overlap_Input.csv", sep="/"), quote=F, row.names=F)

##STEP 6: Create map using all possible markers by combining each species map with the overlap map
#giving the most weight to the overlap map and more weight to the white spruce map due to greater sample size
#load input data (from above)
sitka_lp <- read.csv(file=paste(datadir, "Sitka_LPMerge_Input.csv", sep="/"), header=T)
white_lp <- read.csv(file=paste(datadir, "White_LPMerge_Input.csv", sep="/"), header=T)
combo_lp <- read.csv(file=paste(datadir, "Consensus_Overlap_Input.csv", sep="/"), header=T)

#Need to get a separate dataframe for each linkage group
SitkaLGsall <- split( sitka_lp , f = sitka_lp$LG)
WhiteLGsall <- split( white_lp , f = white_lp$LG)
ComboLGsall <- split( combo_lp , f = combo_lp$LG)

#Then need to remove linkage group column
sitka_input_all <- lapply(SitkaLGsall, function(x) { x["LG"] <- NULL; x })
white_input_all <- lapply(WhiteLGsall, function(x) { x["LG"] <- NULL; x })
combo_input_all <- lapply(ComboLGsall, function(x) { x["LG"] <- NULL; x })

#Need to put linkage groups in the same order for each species
#i.e. if Sitka LG1 is the same as White LG4, the both need to be first in the linkage group list
#this will simplify merging each linkage group
#Aggregate table as key:
#     LG.PS LG.WS   x
# 7     LG4  LG01 240
# 17    LG3  LG02 239
# 26    LG5  LG03 203
# 35    LG9  LG04 198
# 43    LG6  LG05 233
# 56    LG8  LG06 205
# 57    LG1  LG07 324
# 70    LG2  LG08 241
# 78   LG12  LG09 168
# 85   LG11  LG10 141
# 101   LG7  LG11 211
# 104  LG10  LG12 197

#Get white spruce linkage groups in matching order
white_input_all_order <- list(white_input_all[["LG07"]], white_input_all[["LG08"]], white_input_all[["LG02"]], white_input_all[["LG01"]], white_input_all[["LG03"]], white_input_all[["LG05"]], white_input_all[["LG11"]], white_input_all[["LG06"]], white_input_all[["LG04"]], white_input_all[["LG12"]], white_input_all[["LG10"]], white_input_all[["LG09"]])

#AND because R puts numbered linkage groups in numerical order (reading 1:9 incorrectly because there is no zero)
#Fix sitka asn overlap input order (overlap map using Sitka linkage group names)
sitka_input_all_order <- list(sitka_input_all[["LG1"]], sitka_input_all[["LG2"]], sitka_input_all[["LG3"]], sitka_input_all[["LG4"]], sitka_input_all[["LG5"]], sitka_input_all[["LG6"]], sitka_input_all[["LG7"]], sitka_input_all[["LG8"]], sitka_input_all[["LG9"]], sitka_input_all[["LG10"]], sitka_input_all[["LG11"]], sitka_input_all[["LG12"]])
combo_input_all_order <- list(combo_input_all[["LG1"]], combo_input_all[["LG2"]], combo_input_all[["LG3"]], combo_input_all[["LG4"]], combo_input_all[["LG5"]], combo_input_all[["LG6"]], combo_input_all[["LG7"]], combo_input_all[["LG8"]], combo_input_all[["LG9"]], combo_input_all[["LG10"]], combo_input_all[["LG11"]], combo_input_all[["LG12"]])

##Data check
#Check that the Sitka input has only unique markers
#(proving that renaming SNPs and removing duplicate SNPs on GCATs worked correctly)
rows <- list()
unique <- list()
for(x in 1:12){
  rows[[x]] <- data.frame("LG"=paste(x), "Num"=nrow(sitka_input_all[[x]]))
  unique[[x]] <- data.frame("LG"=paste(x), "Num"=length(unique(sitka_input_all[[x]][,1])))
}
rowtot <- do.call(rbind, rows)
unitot <- do.call(rbind, unique)
check_marks <- cbind(rowtot, unitot)
colnames(check_marks) <- c("LG", "Num.Row", "LG", "Num.Uni")
length(check_marks$Num.Row==check_marks$Num.Uni)
#Should be 12
sum(unitot$Num)
#20983

#Check that the ordered LG lists are in the correct order
#The total number of overlapping markers when this is checked by linkage group
#should be equal to the total number of overlapping markers in the two input files (2327)
#check
chk <- list()
for(g in 1:12){
  chk[[g]] <- data.frame("LG"=paste(g), "Overlap"=nrow(sitka_input_all_order[[g]][which(sitka_input_all_order[[g]][,1] %in% white_input_all_order[[g]][,1]),]), "Overlap_Consensus_PS"=nrow(sitka_input_all_order[[g]][which(sitka_input_all_order[[g]][,1] %in% combo_input_all_order[[g]][,1]),]), "Overlap_Consensus_PG"=nrow(white_input_all_order[[g]][which(white_input_all_order[[g]][,1] %in% combo_input_all_order[[g]][,1]),]))
}
chktot <- do.call(rbind, chk)
sum(chktot$Overlap)
#Should be 2327
sum(chktot$Overlap_Consensus_PS)
#Should be 2327
sum(chktot$Overlap_Consensus_PG)
#Should be 2327

#Run LPMerge (as above)
sink(file=paste(outputdir, "LPMerge_Output_wconsensus_weighted.txt", sep="/"))
outlistw_all <- list()
for(x in 1:12){
  print(paste("------------------------","CHR", x, "----------------------------", sep=""))
  nam <- paste("CHR", x, sep="")
  outlistw_all[[nam]] <- LPmerge(list(sitka_input_all_order[[x]], white_input_all_order[[x]], combo_input_all_order[[x]]), max.interval=1:10, weights=c(1, 2, 3))
}
sink()

maxlen <- list()
for(c in 1:12){
  maxlen[[c]] <- data.frame("CHR"=paste("CHR", c, sep=""), "MAXLEN_Sitka"=max(sitka_input_all_order[[c]][,2]),  "MAXLEN_White"=max(white_input_all_order[[c]][,2]), "Average"=mean(c(max(sitka_input_all_order[[c]][,2]), max(white_input_all_order[[c]][,2]))))
}

maxlength <- do.call("rbind", maxlen)

#visualize output of each chromosome
outw <- outlistw_all
chrsw <- list()
for(n in 1:12){
  for(x in 1:10){
    outw[[n]][[x]]$Interval <- rep(x, nrow(outw[[n]][[x]]))
  }
  chrsw[[n]] <- do.call("rbind", outw[[n]])
  chrsw[[n]]$Interval <- as.factor(chrsw[[n]]$Interval)
  print(ggplot(chrsw[[n]], aes(x=consensus, y=Interval))+geom_point()+labs(title=paste("CHR", n))+geom_vline(xintercept=c(maxlength[n,4], max(chrsw[[n]][,3], na.rm=T), max(chrsw[[n]][,4], na.rm=T), max(chrsw[[n]][,5], na.rm=T)), color=c("black", "red", "blue", "green")))
  ggsave(paste(paste(outputdir, "/weighted_plots/", sep=""), "CHR", n, "_allmarks_weighted_withconsensus.png", sep=""), height = 10, width=20, units="cm")
  chrsw[[n]]$LG <- paste("LG", n, sep="")
}

#save output for all intervals for each chromosome together in a single file
chrsw_all <- do.call("rbind", chrsw)
write.csv(chrsw_all, paste(outputdir, "all_markers_output_weighted_withconsensus.csv", sep="/"), row.names = F, quote=F)

##select best interval and remove gaps

#Select best interval
#Determine lengths of the linkage groups in each species
# CHR MAXLEN_Sitka MAXLEN_White  Average
# 1   CHR1      212.424       196.12 204.2720
# 2   CHR2      183.039       161.77 172.4045
# 3   CHR3      198.488       178.94 188.7140
# 4   CHR4      193.762       170.00 181.8810
# 5   CHR5      164.071       157.34 160.7055
# 6   CHR6      172.375       151.28 161.8275
# 7   CHR7      199.389       149.49 174.4395
# 8   CHR8      197.546       164.35 180.9480
# 9   CHR9      163.745       141.97 152.8575
# 10 CHR10      154.818       128.55 141.6840
# 11 CHR11      128.014       124.35 126.1820
# 12 CHR12      142.376       143.89 143.1330

#select best intervals (based on length, RMSE, or both)
#indicate whether graph created above shows gaps in the chromosome
#i.e. a large space with no markers (created due to differences in chromosome length)
##CHR1=2 (RMSE) - gap
##CHR2=2 (RMSE) - gap
##CHR3=7 (RMSE) - gap
##CHR4=2 (RMSE) - gap
##CHR5=2 (RMSE) - no gaps
##CHR6=2 (RMSE) - gaps
##CHR7=2 (RMSE) - gaps
##CHR8=2 (RMSE) - gaps
##CHR9=9 (RMSE) - gaps
##CHR10=2 (RMSE) - gaps
##CHR11=1 (RMSE) - gaps
##CHR12=2 (RMSE) - no gaps

#bring back in file with all intervals
all <- read.csv(paste(outputdir, "all_markers_output_weighted_withconsensus.csv", sep="/"))
#select best intervals from above and save 
best <- all[which(all$LG=="LG1"&all$Interval=="2"|all$LG=="LG2"&all$Interval=="2"|all$LG=="LG3"&all$Interval=="7"|all$LG=="LG4"&all$Interval=="2"|all$LG=="LG5"&all$Interval=="2"
                  |all$LG=="LG6"&all$Interval=="2"|all$LG=="LG7"&all$Interval=="2"|all$LG=="LG8"&all$Interval=="2"
            |all$LG=="LG9"&all$Interval=="9"|all$LG=="LG10"&all$Interval=="2"|all$LG=="LG11"&all$Interval=="1"|all$LG=="LG12"&all$Interval=="2"),]
write.csv(best, file=paste(outputdir, "all_markers_output_weighted_withconsensus_bestchrs.csv", sep="/"))

##STEP 7: Manually remove gaps in linkage groups in excel using graphs from above for loop to inform removal
#NOTE: gaps appear to be formed when markers at terminal end are only on Sitka spruce
#because Sitka has such longer chromosomes

#bring back in map with gaps removed
nogaps <- read.csv(paste(outputdir, "all_markers_output_weighted_withconsensus_bestchrs_nogaps.csv", sep="/"))
#27,052 rows

#visual check
nogap_chk <- ggplot(nogaps, aes(x=consensus, y=LG, color=LG))+geom_point()

#no gaps at end, but note that CHR7 still has gaps within chromosome

##STEP 8A: Calculate map statistics 
chrlen <- aggregate(nogaps$consensus, by=list(nogaps$LG), max)
write.csv(chrlen, paste(outputdir, "Integrated_Map_weighted_withconsensus_NoGaps_ChromosomeLengths.csv", sep="/"))
totallength <- sum(chrlen$x)
#1860.241

#calculate gaps between markers in cM
require(dplyr)
map_order <- nogaps[order(nogaps$LG, nogaps$consensus),]
gaps <- map_order %>%
  group_by(LG) %>%
  mutate(Gaps = consensus - lag(consensus, default = consensus[1]))

gapdf <- as.data.frame(gaps)
gapavg <- aggregate(gapdf$Gaps, by=list(gapdf$LG), mean)
mean(gapavg$x)
#0.06909465
mean(gapdf$Gaps)
#0.06876538
max(gapdf$Gaps)
#16.64
write.csv(gapdf, paste(outputdir, "Integrated_Mapweighted_withconsensus_NoGaps_Gaps.csv", sep="/"))

##STEP 8B: Compare to Sitka and White spruce maps
#Load 3 maps 
#NOTE: Use Sitka LPMerge Input rather than map due to renaming SNPIDs as GCAT
sitka_input <- read.csv(paste(datadir, "Sitka_LPMerge_Input.csv", sep="/"))
white <- read.csv(paste(datadir, "White_Spruce.csv", sep="/"), header=T)
int <- read.csv(paste(outputdir,"all_markers_output_weighted_withconsensus_bestchrs_nogaps.csv", sep="/"), header=T)
int <- int[,-1] #row numbers were saved

#Rename columns for merging
colnames(int)[7] <- "LG.Int"
colnames(int)[2] <- "Position.Int"
colnames(int)[1] <- "SNPID"
colnames(sitka_input)[1:3] <- paste(colnames(sitka_input)[1:3], "PSIn", sep=".")

#check against Sitka map and Sitka input (i.e. one SNP per GCAT)
intin <- merge(sitka_input, int, by.x="Marker.PSIn", by.y="SNPID")#20,857
intin_lg <- aggregate(intin$Position.Int, by = intin[c('LG.PSIn', 'LG.Int')], length)
#all group the same
#check correlation in LGs that match between Sitka input and integrated map
INTINLGs <- split(intin, f = intin$LG.PSIn)
INTINcorr <- list()
for(x in 1:12){
  nam <- paste("LG", x, sep="")
  INTINcorr[[nam]]  <- data.frame("Method"=cor.test(INTINLGs[[x]]$Position.PSIn, INTINLGs[[x]]$Position.Int, method="kendall", ci=FALSE)$method, "Estimate"=cor.test(INTINLGs[[x]]$Position.PSIn, INTINLGs[[x]]$Position.Int, method="kendall", ci=FALSE)$estimate, "Pvalue"=as.numeric(cor.test(INTINLGs[[x]]$Position.PSIn, INTINLGs[[x]]$Position.Int, method="kendall", ci=FALSE)$p.value))
}
INTINcorrs <- do.call("rbind", INTINcorr)
write.csv(INTINcorrs, paste(outputdir, "Integrated_withconsensus_nogaps_SitkaInput_Correlations.csv", sep="/"))
mean(abs(INTINcorrs$Estimate))
# 0.9610085
min(abs(INTINcorrs$Estimate))
# 0.8984562
max(abs(INTINcorrs$Estimate))
#0.9952256

#How many markers added?
nrow(int[-which(int$SNPID %in% sitka_input$Marker.PSIn),])
#6195

#Compare to white Spruce map to integrated map 
intwhite <- merge(white, int, by.x="GCAT", by.y="SNPID")#8,550
intwhite_lg <- aggregate(intwhite$Position.Int, by = intwhite[c('LG.WS', 'LG.Int')], length)
#all markers group the same
#check correlation in LGs that match
INTWSLGs <- split(intwhite, f = intwhite$LG.Int )
INTWScorr <- list()
for(x in 1:12){
  nam <- paste("LG", x, sep="")
  INTWScorr[[nam]]  <- data.frame("Method"=cor.test(INTWSLGs[[x]]$Position.WS, INTWSLGs[[x]]$Position.Int, method="kendall", ci=FALSE)$method, "Estimate"=cor.test(INTWSLGs[[x]]$Position.WS, INTWSLGs[[x]]$Position.Int, method="kendall", ci=FALSE)$estimate, "Pvalue"=as.numeric(cor.test(INTWSLGs[[x]]$Position.WS, INTWSLGs[[x]]$Position.Int, method="kendall", ci=FALSE)$p.value))
}
INTWScorrs <- do.call("rbind", INTWScorr)
write.csv(INTWScorrs, paste(outputdir, "Integrated_withconsensus_nogaps_WhiteSpruce_Correlations.csv", sep="/"))
mean(abs(INTWScorrs$Estimate))
#0.9911118
min(abs(INTWScorrs$Estimate))
#0.9827445
max(abs(INTWScorrs$Estimate))
#0.9992189

#How many markers added?
newwhite <- int[-which(int$SNPID %in% white$GCAT),]
#18519
newGCAT <- newwhite[-which(newwhite$SNPID %in% newwhite$SNPID[grepl("SNP", newwhite$SNPID)]),]
##2798


#Visualize Comparisons
library(ggplot2)

INTPSplot <- ggplot(intin, aes(x=Position.PSIn, y=Position.Int))+geom_point(size=1)+
  ylab("Position in Integrated Map")+ xlab("Position in Sitka Spruce Map")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_grid(LG.Int~LG.PSIn)

ggsave(paste(outputdir, "SitkaMap_IntegratedMap_withconsensus_nogaps.png", sep="/"), INTPSplot, height=6, width=6, units="in")

INTWSplot <- ggplot(intwhite, aes(x=Position.WS, y=Position.Int))+geom_point(size=1)+
  ylab("Position in Integrated Map")+ xlab("Position in White Spruce Map")+
  theme_bw()+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank())+
  theme(panel.grid.major=element_line(colour="white"),panel.grid.minor=element_line(colour="white"))+
  theme(panel.spacing.x=unit(0.01, "lines"),panel.spacing.y=unit(0.01, "lines"))+
  facet_grid(LG.Int~LG.WS)

ggsave(paste(outputdir, "WhiteMap_IntegratedMap_withconsensus_nogaps.png", sep="/"), INTWSplot, height=6, width=6, units="in")

#Overlay plots
colnames(sitka_input)[1] <-"LG"
colnames(sitka_input)[3] <- "Position"
colnames(int)[2] <- "Position"
colnames(int)[7] <- "LG"
colnames(white)[4] <- "LG"
colnames(white)[5] <- "Position"
white[white == 'LG01'] = 'LG1'
white[white == 'LG02'] = 'LG2'
white[white == 'LG03'] = 'LG3'
white[white == 'LG04'] = 'LG4'
white[white == 'LG05'] = 'LG5'
white[white == 'LG06'] = 'LG6'
white[white == 'LG07'] = 'LG7'
white[white == 'LG08'] = 'LG8'
white[white == 'LG09'] = 'LG9'
sitka_input$LG <- factor(sitka_input$LG, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
int$LG <- factor(int$LG, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
white$LG <- factor(white$LG, levels=c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10", "LG11", "LG12"))
overlay <- ggplot(sitka_input, aes(x=Position, y=LG))+geom_point(color="#0066CC")+
  geom_point(data=white, color="#339966")+geom_point(data=int, color="#FFCC33")

ggsave(paste(outputdir, "AllMaps_withconsensus_nogaps.png", sep="/"), overlay, height=6, width=8, units="in")
