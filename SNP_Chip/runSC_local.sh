#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=SepChrom
#SBATCH --time=72:00:00


YOURDIR=#DataFiles
cd $YOURDIR

LEPMAPDIR=#LepMap3Directory/bin
DAT=CombQC_NRf.call
mkdir SC_OUT
OUTDIR=#Output/SC_OUT



for j in $(seq 20 90)
do
java -cp $LEPMAPDIR SeparateChromosomes2 data=$DAT lodLimit=$j sizeLimit=100 families=LM2 > $OUTDIR/SepChrom_Lod$j.txt
done
