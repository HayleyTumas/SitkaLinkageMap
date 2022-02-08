#Need to create fasta file for all mapped RADSeq sequences to blast against white spruce transcriptome

install.packages("seqinr") #v4.2-5
library(seqinr)

datafiles <- #Path to DataFiles
output <- #Path to output

#STEP1 1: Load RADSeq sequence data
#And data on which SNPs were identified in 3 families from the Sitka spruce breeidng population
#Note: Only families 1 and 2 were used to construct the linkag map
##Load RADSeq data from Joanna Ilska
#NOTE: Will need to unzip tags file using gunzip
tag <- read.table(paste(datafiles, "catalog.tags.tsv", sep=""), sep='\t', header=F, skip=1)
f1 <- read.table(paste(datafiles, "FAM1_SNPs.txt", sep=""), header=F)
f2 <- read.table(paste(datafiles, "FAM2_SNPs.txt", sep=""), header=F)
f3 <- read.table(paste(datafiles, "FAM3_SNPs.txt", sep=""), header=F)

#STEP 2: Get list of SNPs identified in the 3 families
#Note that "Position" is not a unique SNP identifier
#And that there are up to 3 SNPs per locus
map_snp <- rbind(f1, f2, f3)
#get rid of duplicated
map_snp_uni <- map_snp[-which(duplicated(map_snp$V3)),]

#STEP 3: Subset the tag file
#V2 in the tag file is the locus ID or locus (POS in the vcf)
#V2 in the snp file is also the locus ID, multiple SNPs as above. 
#V3 in the snp file is the position of the SNP on the locus sequences (1-48)
#the pos/chr has been swapped in the vcf
#so locus ID=POS and position on locus=CHR
#makes sense as more duplicates of CHR and only goes upt o 48
#The fasta file will uses the POS from the vcf file
#but keep the link file from the vcf file too to link the SNPID in the map results
tag_mapped <- tag[which(tag$V2 %in% map_snp_uni$V2),]

#STEP 4: Subset this to include only SNPs on 'final' map
#And create a fasta file
fullmap <- read.csv(paste(datafiles, "RADSNP_SCLM2_Cook4NARC_Map_SNPID.csv", sep=""))
radmap <- fullmap[which(fullmap$Source=="RAD"),] #NOTE: All RADSeq SNPs contain "SNP" in ID (which Chip SNPs do not) - can retrieve via grepl
#check the number of positions (loci) of the RADSeq sequences on the map
length(unique(radmap$POS))
#16975 - which matches the total number of RADSeq SNPs (nrow(radmap))
#which means that only unique loci are on the map (as was intended)
#Subset all SNPs 
tag_mapped <- tag[which(tag$V2 %in% radmap$POS),] #16975, good 
tag_to_fasta <- tag_mapped[,c("V2", "V6")]
#since we only have one SNP per locus, it would be more informative to have the SNPID in the fasta
tagID <- merge(tag_to_fasta, radmap[,c("SNPID","POS")], by.x="V2", by.y="POS")
#reconfigure so we've got a good ID
tagID$ID <- paste(tagID$SNPID, (paste("CHR",tagID$V2, sep="_")), sep=".")
tagform <- data.frame("ID"=tagID$ID, "seq"=tagID$V6)
fa = character(2 * nrow(tagform))
fa[c(TRUE, FALSE)] = sprintf("> %s", tagform$ID)
fa[c(FALSE, TRUE)] = as.character(tagform$seq)
writeLines(fa, paste(output "RADSNP_mapped_radseq.fasta", sep=""))
