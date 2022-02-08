#Subset Norway spruce genome to only mapped SNPs
#Use the Picea abies genome 1.0
#downloaded from the FTP on this site:http://congenie.org/sequence_search

datafiles <- #Path to DataFiles/

library(seqinr)

#load file of all mapped SNPs from Bernhardsson et al. 2019
df <- read.csv(paste(datafiles, "Bernhardssonetal2019.csv", sep=""))

#STEP 1: Create file of sequences to extract in seqtk
#create text file of seq names
#multiple SNPs on one contig - only want unique contigs (Scaffolds)
unidf <- df[-which(duplicated(df$Scaffold)),]
length(unique(df$Scaffold))
uniscaff <- as.data.frame(unidf$Scaffold)
write.table(uniscaff, file=paste(datafiles, "PA_SNPs.txt", sep=""), sep = "\t", col.names = F, row.names = F, quote = F)
#subset using seqtk (see Data_Prep.sh)
#Double check fasta file
pasub <- read.fasta(paste(datafiles, "PA_SNP.fasta", sep=""))


#let's look at the blast results 
#used e-50 on forward, e-100 on reverse because why make it easy
colnames <- c("qseqid", "sseqid", "sacc", "qlen", "slen", "nident", "length", "mismatch", "gaps", "bitscore", "evalue")
gfor <- read.csv(paste(datafiles, "GCAT_Pabies.csv", sep=""))
grev <- read.csv(paste(datafiles, "GCAT_Pabies_REV.csv", sep=""))
pfor <- read.csv(paste(datafiles, "PSIR_Pabies.csv", sep=""))
prev <- read.csv(paste(datafiles, "PSIR_Pabies_REV.csv", sep=""))
colnames(gfor) <- colnames
colnames(grev) <- colnames
colnames(pfor) <- colnames
colnames(prev) <- colnames

#get only matches <e-100
gfor100 <- gfor[which(gfor$evalue<=1E-100),]
pfor100 <- pfor[which(pfor$evalue<=1E-100),]

#now get only those matches that also have a reverse
#first get only unique matches
gfor_order <- gfor100[order(gfor100$qseqid, gfor100$sseqid),]
gforu <- gfor_order[!duplicated(gfor_order[,c('qseqid', 'sseqid')]),]

pfor_order <- pfor100[order(pfor100$qseqid, pfor100$sseqid),]
pforu <- pfor_order[!duplicated(pfor_order[,c('qseqid', 'sseqid')]),]

grev_order <- grev[order(grev$qseqid, grev$sseqid),]
grevu <- grev_order[!duplicated(grev_order[,c('qseqid', 'sseqid')]),]

prev_order <- prev[order(prev$qseqid, prev$sseqid),]
prevu <- prev_order[!duplicated(prev_order[,c('qseqid', 'sseqid')]),]

#now get only those with matches in both
gforu$ID <- paste(gforu$qseqid, gforu$sseqid)
grevu$ID <- paste(grevu$sseqid, grevu$qseqid)
gmatch <- merge(gforu, grevu, by="ID")

pforu$ID <- paste(pforu$qseqid, pforu$sseqid)
prevu$ID <- paste(prevu$sseqid, prevu$qseqid)
pmatch <- merge(pforu, prevu, by="ID")

write.csv(gmatch, file=paste(datafiles, "GCAT_PAbies_e-100_revfor.csv", sep=""))
write.csv(pmatch, file=paste(datafiles, "PSIR_PAbies_e-100_revfor.csv", sep=""))
