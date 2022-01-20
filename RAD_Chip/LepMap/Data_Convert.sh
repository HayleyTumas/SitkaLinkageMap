
#On local machine for first attempt with GenomeStudio data
#activate environment with Java
conda activate LepMap

#Lepmap directory
LMdir=#LepMap3directory/bin

#go to directory with files
cd #DataFiles

#Load and convert  data filtered for missingness (ind=0.15, SNP=0.2) and maf (0.2) with no zeroed SNPs or positive controls (5533 markers)
java -cp $LMdir ParentCall2 data=RAD_SNPInds_ped.txt vcfFile=RAD_SNP_2FAMcomb.vcf removeNonInformative=1 >RAD_SNP.call
#gives 24702 informative markers

#filter data
java -cp $LMdir Filtering2 data=RAD_SNP.call dataTolerance=0.01 > RAD_SNPf.call
#24702

#now by ByFamily
java -cp $LMdir ParentCall2 data=RAD_SNP_LM1ped.txt vcfFile=RAD_SNP_LM1.recode.vcf removeNonInformative=1 >RAD_SNP_LM1.call

#filter data
java -cp $LMdir Filtering2 data=RAD_SNP_LM1.call dataTolerance=0.01 > RAD_SNP_LM1f.call

#FAM2
java -cp $LMdir ParentCall2 data=RAD_SNP_LM2ped.txt vcfFile=RAD_SNP_LM2.recode.vcf removeNonInformative=1 >RAD_SNP_LM2.call

#filter data
java -cp $LMdir Filtering2 data=RAD_SNP_LM2.call dataTolerance=0.01 > RAD_SNP_LM2f.call

#And now with poor aligning SNPs filtered out
#using only mapped markers minus markers greater than 4/n
java -cp $LMdir ParentCall2 data=RAD_SNPInds_ped.txt vcfFile=RAD_SNP_2FAM_Cook4N.recode.vcf removeNonInformative=1 >RAD_SNP_RPCook4N.call
#gives 21571 informative markers - no markers lost, thats a good sign

#filter data
java -cp $LMdir Filtering2 data=RAD_SNP_RPCook4N.call dataTolerance=0.01 > RAD_SNP_RPCook4N_f.call
