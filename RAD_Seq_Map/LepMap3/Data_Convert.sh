
#On local machine for first attempt with GenomeStudio data
#activate environment with Java
conda activate LepMap

#Lepmap directory
LMdir=#LepMap3Directory

#go to directory with files
YOURDIR=#DataFiles
cd $YOURDIR


java -cp $LMdir ParentCall2 data=BCF_pedtR.txt vcfFile=RAD_comb2_QC_BestSNP.recode.vcf removeNonInformative=1 >RADQC_2FAMSS.call
#19529 informative markers

#filter data
java -cp $LMdir Filtering2 data=RADQC_2FAMSS.call dataTolerance=0.01 > RADQC_2FAMSSf.call

#By Family
java -cp $LMdir ParentCall2 data=RADQC_FAM1_pedtR.txt vcfFile=fam1_qcedSNPs.recode.vcf removeNonInformative=1 >RADQC_FAM1.call
