#!/bin/bash

YOURDIR=#Output

cd $YOURDIR

mkdir Top_OM

echo "Chromosome Likelihoods" > mapinfo.txt
echo 'C=Chromosome, I=Iteration, N.SNP=Number SNPs, LIKE=Likelihood, LEN=Map Length' >meta.txt
cat meta.txt >> mapinfo.txt
rm meta.txt


for s in $(seq 1 12)
do
  cd $YOURDIR/OM_CHR${s}
  echo 'C'" "'I'" "'N.SNP'" "'LIKE'" "" ""   "" "" "'LENDAD'" "" "'LENMOM' > mapinfo_CH${s}.txt
  MAX=-1000000000
  BEST=1
  for j in $(seq 1 5)
  do
    nsnp=$(cat CHR${s}_ITER${j}.txt | awk 'END{print NR-4}')
    echo ${s}" "${j}" "${nsnp} > out.chr
    cat CHR${s}_ITER${j}.txt | grep -o 'likelihood = .*'| awk '{print $3}'> out.likl
    cat CHR${s}_ITER${j}.txt | tail -1 | awk 'NR==1{print $2,$3}' > out.len
    paste -d'  ' out.chr out.likl out.len > out.all
    cat out.all >> mapinfo_CH${s}.txt
    LL=$(cat out.likl)
    if [ 1 -eq "$(echo "${LL} > ${MAX}" | bc)" ]; then MAX=${LL}; BEST=${j}; fi
    rm out.*
  done
  cp CHR${s}_ITER${BEST}.txt $YOURDIR/Top_OM
  echo '----------Chromosome '${s}' BEST='$BEST'-------' > Top.txt
  cat Top.txt mapinfo_CH${s}.txt >> $YOURDIR/mapinfo.txt
  rm Top.txt
  BEST=1
  rm ./-1000000000
done
