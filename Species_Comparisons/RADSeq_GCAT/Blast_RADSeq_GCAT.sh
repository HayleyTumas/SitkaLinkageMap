YOURDIR=#DataFiles
cd $YOURDIR
#GCAT transcriptome fasta file from Rigault et al. 2011
#RADSeq fasta file created with Data_Prep R code

#open blast environment containing blast-plus v2.2.31
conda activate BLAST

#created new file of RADseq loci mapped in 'final' RADSNP map
#do same thing, but I also want pident to compare to the norway spruce paper stats
#create databases with full transcriptome and gene catalog
makeblastdb -in RADSNP_mapped_radseq.fasta -dbtype nucl -out mappedRAD_db
makeblastdb -in GCAT_WS-3.3.cluseq.fa -dbtype nucl -out GCAT_Full_db

#Full Gene Catalog
blastn -query RADSNP_mapped_radseq.fasta -db GCAT_Full_db -evalue 1e-10 -outfmt "10 qseqid sseqid qlen slen pident nident length mismatch gaps bitscore evalue" -out RADSeq_FullGCAT.csv
blastn -query GCAT_WS-3.3.cluseq.fa -db mappedRAD_db -evalue 1e-10 -outfmt "10 qseqid sseqid qlen slen pident nident length mismatch gaps bitscore evalue" -out RADSeq_FullGCAT_REV.csv
