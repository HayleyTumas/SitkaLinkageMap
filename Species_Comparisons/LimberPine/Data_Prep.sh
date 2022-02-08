YOURDIR=/mnt/c/Users/dops0764/Documents/BLAST/Conifers/Pinusflexilis

cd $YOURDIR

#activate environment containing seqtk
conda activate FASTAsub

#Transcriptome fasta file was sent by Liu et al. 2019 corresponding author Jun-Jun Liu
#check to see if all mapped SNPs are included in this fasta file
grep -e ">" Mapped-g9,520-cds-nt.fasta >namessent.txt

#Analyzed list in R - compared to map file from publication (S5 table)
#Found 100 missing based on SNP list in table S5 and 6 in fasta that are not in table
#subset to exclude 6 not in table
seqtk subseq Mapped-g9,520-cds-nt.fasta PF_SNPs_File.txt > PF_Mapped_SNP.fasta
