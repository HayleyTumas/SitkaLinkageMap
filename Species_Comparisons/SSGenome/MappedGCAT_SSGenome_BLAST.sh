YOURDIR=#DataFiles
cd $YOURDIR
OUTDIR=#Ouptut
#created fasta of GCAT matches with RAD and SNP chip that were mapped
#found in R code for analyzing RADSeq to GCAT Blast results

#open blast environment containing blast-plus v2.2.31
conda activate BLAST

#create database with mapped GCAT
makeblastdb -in GCAT_Mapped_RADChip.fasta -dbtype nucl -out mappedGCAT_db
#create database for Sitka spruce Sitka spruce genome
#NOTE: Genome file not found here, too large
#Genome can be found on GenBank: https://www.ncbi.nlm.nih.gov/assembly/GCA_010110895.1/
#Assembly accession #:GCA_010110895.1
#Uploaded by BC Cancer Agency 2020/01/09
makeblastdb -in GCA_010110895.1_Q903_v1_genomic.fna -dbtype nucl -out ssgenome_db

blastn -query GCAT_Mapped_RADChip.fasta -db ssgenome_db -num_threads 16 -evalue 1e-10 -outfmt "10 qseqid sseqid sacc qlen slen pident nident length mismatch gaps bitscore evalue" -out $OUTDIR/RADSeq_SSGenome.csv
blastn -query GCA_010110895.1_Q903_v1_genomic.fna -db mappedGCAT_db -num_threads 16 -evalue 1e-10 -outfmt "10 qseqid sseqid sacc qlen slen pident nident length mismatch gaps bitscore evalue" -out $OUTDIR/RADSeq_SSGenome_REV.csv
#note: moved to SSgenome directory to avoid copying genome file again
