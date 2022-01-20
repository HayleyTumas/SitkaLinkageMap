#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=SepChrom
#SBATCH --time=72:00:00

module load java

YOURDIR=#Output
DATDIR=#DataFiles
LEPMAPDIR=#LepMap3directory/bin
cd $YOURDIR

DAT=$DATDIR/RAD_SNP_LM2f.call
mkdir SC_OUT
OUTDIR=$YOURDIR/SC_OUT



for j in $(seq 15 90)
do
java -cp $LEPMAPDIR SeparateChromosomes2 data=$DAT lodLimit=$j sizeLimit=100 > $OUTDIR/SepChrom_Lod$j.txt
done
