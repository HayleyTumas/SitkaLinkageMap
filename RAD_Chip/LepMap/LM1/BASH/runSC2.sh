#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=SepChrom
#SBATCH --time=72:00:00

module load java

YOURDIR=#Output
DATDIR=#DataFiles
LEPMAPDIR=#LepMap3directory/bin
cd $YOURDIR

DAT=$DATDIR/RAD_SNP_LM1f.call
mkdir SC_OUT2
OUTDIR=$YOURDIR/SC_OUT2
map42=SepChrom_Lod42.txt

cp $YOURDIR/SC_OUT/$map42 $YOURDIR

for j in $(seq 25 55)
do
java -cp $LEPMAPDIR SeparateChromosomes2 data=$DAT lodLimit=$j sizeLimit=100 map=$map42> $OUTDIR/SepChrom_Lod$j.txt
done
