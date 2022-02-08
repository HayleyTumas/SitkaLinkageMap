#Create fasta file of mapped sequences for limber pine

datafiles <- #Path to DataFiles

library(seqinr)

#load file of all mapped SNPs from Liu et al 2019 (doi: 10.1111/tpj.14270)
df <- read.csv(paste(datafiles, "pflexilis_map.csv", sep=""))

#Checked the Picea abies genome 1.0 available on congenie:http://congenie.org/sequence_search
#Did not have all of the scaaffolds represented 
#Contacted corresponding author Jun-Jun Liu
#received Mapped-g9,520-csd-nt.fasta file of mapped sequences

##STEP 1: create file of mapped sequence names to check they match fasta file
#Multiple SNPs per contig (scaffold or sequence) - find unique contigs
unidf <- df[-which(duplicated(df$Transcript.ID))]
length(unique(df$Transcript.ID))
uniscaff <- as.data.frame(df$Transcript.ID)
write.table(uniscaff, file="PF_SNPs.txt", sep = "\t", col.names = F, row.names = F, quote = F)

##STEP 2: Compare contig names in publication map file (table S5) to fasta sent by Liu
#fasta sequence names were extracted using seqtk (check accompanying script)
#load sequence names extracted by seqtk
fastanames <- read.table("namessent.txt")
#reformat to match format in SNP table
fastanames$ID <-gsub("[\\>]","",fastanames[,1])
miss <- as.data.frame(uniscaff[-which(uniscaff$`df$Transcript.ID` %in% fastanames$ID),])
#100 sequences are in the publication map file (table S5), but not in fasta file sent by Liu
missing <- as.data.frame(fastanames[-which(fastanames$ID %in% uniscaff$`df$Transcript.ID`),])
#6 sequences are in the fasta file sent by Liu, but not in the publication map file (table S5)
#Want to subset fasta file to only include sequences that appear in map file
#(blast matches to sequences not mapped would be useless information)
#create list of SNPs in both the fasta file and map file
match <- as.data.frame(fastanames[which(fastanames$ID %in% uniscaff$`df$Transcript.ID`),])
matches <- as.data.frame(match$ID)
write.table(matches, file="PF_SNPs_File.txt", sep = "\t", col.names = F, row.names = F, quote = F) 

##STEP 3: Subset fasta file to include only mapped sequences in publication map file
#use seqtk - see accompanying code - output: PF_Mapped_SNP.fasta
#Double check number of sequences here
pf <- read.fasta(paste(datafiles, "PF_Mapped_SNP.fasta", sep=""))
