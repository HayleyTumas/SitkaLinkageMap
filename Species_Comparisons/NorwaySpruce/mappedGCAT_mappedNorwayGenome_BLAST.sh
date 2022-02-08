YOURDIR=$DataFiles
cd $YOURDIR
#created fasta of GCAT matches with RAD and SNP chip that were mapped

#open blast environment containing blast-plus v2.2.31
conda activate BLAST

#create databases with mapped GCAT
makeblastdb -in GCAT_Mapped_RADChip.fasta -dbtype nucl -out mappedGCAT_db
makeblastdb -in PA_SNP.fasta -dbtype nucl -out mappedPA_db
#these files represent GCAT sequences that were mapped with either RAD or Chip SNPs
#and picea albies sequences that were mapped
#for the PA_SNP file check the NorwaySpruce folder in Literature-Other Maps-Spruce
#the code for the GCAT file is in the RAD seq R code in the RADSeq folder

blastn -query GCAT_Mapped_RADChip.fasta -db mappedPA_db -num_threads 16 -evalue 1e-10 -outfmt "10 qseqid sseqid sacc qlen slen pident nident length mismatch gaps bitscore evalue" -out OurMap_NorwaySpruce.csv
blastn -query PA_SNP.fasta -db mappedGCAT_db -num_threads 16 -evalue 1e-10 -outfmt "10 qseqid sseqid sacc qlen slen pident nident length mismatch gaps bitscore evalue" -out OurMap_NorwaySpruce_REV.csv
