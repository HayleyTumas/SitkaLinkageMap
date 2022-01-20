#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=JoinSingles
#SBATCH --time=05:00:00


BESTLOD=28
YOURDIR=#Output
DATDIR=#DataFiles
LEPMAPDIR=#LepMap3directory/bin
MAP=$YOURDIR/SC_OUT/SepChrom_Lod$BESTLOD.txt
DAT=$DATDIR/RAD_SNP_RPCook4N_f.call

cd $YOURDIR

mkdir JS_OUT
OUTDIR=$YOURDIR/JS_OUT


for j in $(seq 2 43)
do
java -cp $LEPMAPDIR JoinSingles2All map=$MAP data=$DAT lodLimit=$j iterate=1 lodDifference=10> $OUTDIR/JS_Lod$j.txt
done
