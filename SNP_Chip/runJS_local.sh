#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=JoinSingles
#SBATCH --time=05:00:00


BESTLOD=41
YOURDIR=#Output
DATDIR=#DataFiles
cp $YOURDIR/SC_OUT/SepChrom_Lod$BESTLOD.txt $YOURDIR
MAP=SepChrom_Lod$BESTLOD.txt
DAT=$DATDIR/CombQC_NRf.call

cd $YOURDIR

LEPMAPDIR=#LepMap3Directory/bin
mkdir JS_OUT_new
OUTDIR=$YOURDIR/JS_OUT_new


for j in $(seq 5 50)
do
java -cp $LEPMAPDIR JoinSingles2All map=$MAP data=$DAT lodLimit=$j iterate=1 lodDifference=10> $OUTDIR/JS_Lod$j.txt
done
