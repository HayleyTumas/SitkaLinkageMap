#!/bin/bash

YOURDIR=#Output/SC_OUT
cd $YOURDIR

TOTMARK=3720

echo 'The number of markers per Linkage group out of '${TOTMARK}'' > nummarkperLG.txt

for i in $(ls -1 $YOURDIR*)
do
lod=$(grep -o 'lodLimit=.* ' ${i})
echo '<<==== '${lod}' ====>>' > lod.txt
cat ${i} | sort -n | uniq -c > nLG.txt
cat lod.txt nLG.txt > lodnLG.txt
cat lodnLG.txt >> nummarkperLG.txt

rm nLG.txt lod.txt lodnLG.txt

done
