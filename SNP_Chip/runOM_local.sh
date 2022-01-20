#!/bin/bash

#SBATCH --nodes=16
#SBATCH --job-name=OrderMark
#SBATCH --time=120:00:00

module load java

BESTLOD=11
YOURDIR=#Output
DATDIR=#Datafiles
cd $YOURDIR
cp $YOURDIR/JS_OUT_new/JS_Lod$BESTLOD.txt $YOURDIR
MAP=JS_Lod$BESTLOD.txt
DAT=$DATDIR/CombQC_NRf.call



LEPMAPDIR=#LepMap3Directory/bin


for j in $(seq 1 12)
do
  mkdir OM_CHR${j}
  OUTDIR=$YOURDIR/OM_CHR${j}
  for m in $(seq 1 5)
  do
    java -cp $LEPMAPDIR OrderMarkers2 map=$MAP data=$DAT chromosome=${j} useKosambi=1 numThreads=16 outputPhasedData=1 sexAveraged=1>$OUTDIR/CHR${j}_ITER${m}.txt
  done
done
