YOURDIR=#DataFiles
cd $YOURDIR

#open blast environment containing blast-plus v2.2.31
conda activate BLAST

#create databases with mapped GCAT and mapped limber pine
#Code for creating mapped GCAT fasta is in the RADSeq to GCAT folder
#limber pine fasta file in Data_Prep code - original fasta sent by Liu et al. 2019 corresponding author
makeblastdb -in GCAT_Mapped_RADChip.fasta -dbtype nucl -out mappedGCAT_db
makeblastdb -in PF_Mapped_SNP.fasta -dbtype nucl -out mappedPF_db

blastn -query GCAT_Mapped_RADChip.fasta -db mappedPF_db -num_threads 16 -evalue 1e-10 -outfmt "10 qseqid sseqid sacc qlen slen pident nident length mismatch gaps bitscore evalue" -out OurMap_LimberPine.csv
blastn -query PF_Mapped_SNP.fasta -db mappedGCAT_db -num_threads 16 -evalue 1e-10 -outfmt "10 qseqid sseqid sacc qlen slen pident nident length mismatch gaps bitscore evalue" -out OurMap_LimberPine_REV.csv
