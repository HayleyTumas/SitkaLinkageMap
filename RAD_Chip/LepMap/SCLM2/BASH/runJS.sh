#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=JoinSingles
#SBATCH --time=05:00:00

module load java

BESTLOD=30
YOURDIR=#Output
DATDIR=#DataFiles
LEPMAPDIR=#LepMap3directory/bin
MAPDIR=#Output/SC_OUT
cd $MAPDIR
cp SepChrom_Lod$BESTLOD.txt $YOURDIR
MAP=SepChrom_Lod$BESTLOD.txt
DAT=$DATDIR/RAD_SNPf.call

cd $YOURDIR

mkdir JS_OUT
OUTDIR=$YOURDIR/JS_OUT


for j in $(seq 10 45)
do
java -cp $LEPMAPDIR JoinSingles2All map=$MAP data=$DAT lodLimit=$j iterate=1 lodDifference=10> $OUTDIR/JS_Lod$j.txt
done
