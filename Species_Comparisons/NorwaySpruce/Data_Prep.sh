YOURDIR=#DataFiles

cd $YOURDIR

#create new environment for seqtk ()
conda create --name FASTAsub

conda install seqtk -n FASTAsub

conda activate FASTAsub

#NOTE: This is Picea abies genome 1.0
#downloaded from the FTP on this site:http://congenie.org/sequence_search
seqtk subseq Pabies1.0-genome.fa.gz PA_SNPs.txt > PA_SNP.fasta
