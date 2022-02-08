
#activate environment with Java
conda activate LepMap

#Lepmap directory
LMdir=#your directory/LepMap/binary+code/bin

#go to directory with files
YOURDIR=#DataFiles
cd $YOURDIR

#Load and convert  data filtered for missingness (ind=0.15, SNP=0.2) and maf (0.2) with no zeroed SNPs or positive controls (5533 markers)
java -cp $LMdir ParentCall2 data=Comb_ped_NoReps.txt vcfFile=Comb_QC_Noreps.recode.vcf removeNonInformative=1 >CombQC_NR.call
#gives 5194 informative markers

#further filter the data in LepMap3
java -cp $LMdir Filtering2 data=CombQC_NR.call dataTolerance=0.01 > CombQC_NRf.call
#5194
