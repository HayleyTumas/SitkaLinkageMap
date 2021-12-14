#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=SepChrom
#SBATCH --time=72:00:00

module load java

YOURDIR=#Output
cd $YOURDIR
DATDIR=#DataFiles

LEPMAPDIR=#LepMap3Directory/bin
DAT=$DATDIR/RADQC_2FAMSSf.call
mkdir SC_OUT
OUTDIR=$YOURDIR/SC_OUT



for j in $(seq 15 95)
do
java -cp $LEPMAPDIR SeparateChromosomes2 data=$DAT lodLimit=$j sizeLimit=100 families=LM2 > $OUTDIR/SepChrom_Lod$j.txt
done
