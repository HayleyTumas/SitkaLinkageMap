#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=JoinSingles
#SBATCH --time=05:00:00

module load java

YOURDIR=#Output
DATDIR=#DataFiles
LEPMAPDIR=#LepMap3directory/bin
BESTLOD=30
MAPDIR=$YOURDIR/SC_OUT
cd $MAPDIR
cp SepChrom_Lod$BESTLOD.txt $YOURDIR
MAP=SepChrom_Lod$BESTLOD.txt
DAT=RAD_SNP_LM2f.call

cd $YOURDIR

mkdir JS_OUT
OUTDIR=$YOURDIR/JS_OUT


for j in $(seq 5 40)
do
java -cp $LEPMAPDIR JoinSingles2All map=$MAP data=$DAT lodLimit=$j iterate=1 lodDifference=10> $OUTDIR/JS_Lod$j.txt
done
