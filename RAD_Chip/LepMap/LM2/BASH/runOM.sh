#!/bin/bash

#SBATCH --nodes=16
#SBATCH --job-name=OrderMark
#SBATCH --time=120:00:00

module load java

BESTLOD=8
YOURDIR=#Output
DATDIR=#DataFiles
LEPMAPDIR=#LepMap3directory/bin
cd $YOURDIR
MAP=JS_Lod$BESTLOD.txt
DAT=$DATDIR/RAD_SNP_LM2f.call

cp $YOURDIR/JS_OUT/$MAP $YOURDIR



for j in $(seq 1 12)
do
  mkdir OM_CHR${j}
  OUTDIR=$YOURDIR/OM_CHR${j}
  for m in $(seq 1 5)
  do
    java -cp $LEPMAPDIR OrderMarkers2 map=$MAP data=$DAT chromosome=${j} useKosambi=1 numThreads=16 outputPhasedData=1 sexAveraged=1>$OUTDIR/CHR${j}_ITER${m}.txt
  done
done
