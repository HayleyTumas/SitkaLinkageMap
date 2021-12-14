# SitkaLinkageMap
Code and accompanying datasets for analyses in Tumas et al. 202X

NOTE: Coding was either completed in Ubuntu for windows using miniconda through Anaconda3 or in R v4.02 using R Studio. LepMap3 v0.2 (4 Dec 2020) (downloaded from SourceForge: https://sourceforge.net/p/lep-map3/wiki/LM3%20Home/) was used to make all final maps that used family 2 for grouping during Separate Chromosomes. Other previous maps that did not use family 2 for grouping would have used the previous version of LepMap3 that does not have this capability. This shows the steps to achieve the maps that were used in the manuscript. Other maps were made, either before proper QC, before we knew to use LM2 during Separate Chromosomes, or individual family maps in RAD-Seq or Chip datasets. Code is in parentheses and output datafiles or directories are in brackets. 
A.	RAD-Seq data
1.	Data (Filter_032020) handed over by Joanna Ilska Warner (code notebook for steps)[ fam1_qcedSNPs.recode.vcf and fam3_qcedSNPs.recode.vcf].
2.	Merged data files for families 1 and 2 using BCFtools from Samtools (merge_VCF.sh)[ RAD_comb2_QC.vcf].
3.	Needed to get single SNP per 48 bp locus. Used vcftools code to find SNP with highest read depth per locus and used this SNP. Read output file in R to find SNP with highest read depth per locus. Subset VCF keeping only “best” SNP per locus using vcftools (snp_depth.sh) [RAD_comb2_QC_BestSNP.recode.vcf]. Note: POS corresponds to locus and CHR corresponds to position on 48bp locus in RADSeq files. 
4.	Used vcfR in R to get individuals (LM_RadSeq.R) [BCF_vcf_inds.csv] and transposed in R (transpose_code) to create pedigree file [BCF_pedtR.txt]. Note: individuals the same as in RAD_comb2_QC.vcf so used individual file and pedigree file from this file.
5.	Converted to LepMap format and further filtered (Data_Convert.sh) [RADQC_2FAMSSf.call].
6.	Using updated LepMap, ran Separate Chromosomes on the ARC cluster using family 2 for grouping (families=LM2) across lodLimits 15-95 (runSC.sh)[SC_OUT]. Selected visually 20 as the lodLimit that assigned the most markers (0 has lowest number) to 12 linkage groups (chk_chr.sh)[nummakerperLG.txt]. 
7.	Ran Join Singles on ARC cluster (runJS.sh) across lodLimits 10-45 [JS_OUT]. Selected 10 as the lodLimit that assigned the most markers (chk_chr_JS.sh) [nummakerperLG.txt].
8.	Ran Order Markers on ARC cluster (runOM.sh) for 5 iterations across each linkage group. Selected the iteration with the highest log likelihood (checklikelihood_CHR.sh)[Top_OM]. 
9.	Combined all files and merged with file (getsnpscall.sh) that matches SNP order in call file to the locus (called POS in the RADSeq files as the unique identifier) [2FAMSS_callfile.txt] and the file that matches the locus to the SNPID [2FAMSS_SNPIDs.txt] (combine_OM_files.R) [RADQC_2FAMSSSCLM2_Map_SNPID].
B.	SNP Chip data
1.	Converted output excel files from GenomeStudio to vcf files using a combination of R (illumine_to_lgen.R) and Plink (Plink_to_vcf.sh). Note: a unique CHR and position number was assigned so that SNPs could be uniquely identified in LepMap call files. Chromosome is 1 for all SNPs (no 1 chromosome in the RADSeq files) and Position is just the row number. Row number was assigned using a list of ALL SNPs as a key then combined with individual files so that SNPs would have the same chr/pos in the 2 separate family files [Combine_ChrKey: LM1_NoZeroSNP.vcf and LM2_NoZeroSNP.vcf].  
2.	Filtered the data in Plink based on missing data (by locus), missingness (by individual) and minor allele frequency (0.2) (Plink_Filtering) [LM1QC.vcf and LM2QC.vcf].
3.	Combine data using Plink (merge_VCF.sh) [Comb_QC.vcf].
4.	Turns out there are some replicate individuals in Family 1. This is because some DNA samples were believed to have failed and were re-extracted. In some instances, the failed samples were not removed from the plate and they did amplify. Removed these from the combined file using VCFtools (no_reps)[ Comb_QC_noreps.recode.vcf].
5.	Create pedigree file using vcfR to get individual names, adding columns for LepMap3 by hand, and then transposing for LepMap3 in R. 
6.	Run in LepMap3. Convert the data to LepMap3 format and run further filtering (Data_Convert.sh) [CombQC_NRf.call].
7.	Run Separate Chromosomes on the local machine across lodLimits 20-90 and using family 2 for grouping (runSC_local.sh) [SC_OUT]. Determine 41 as best lodLimit with fewest unassigned markers and 12 linkage groups (chk_chr_local.sh).
8.	Run Join Singles on the local machine across lodLimits 5-50 (runJS_local.sh)[JS_OUT]. Determine 11 to be best LodLimit (chk_chr_JS_local.sh).
9.	Run Order Markers across 12 linkage groups with 5 iterations each (runOM_local.sh). Determine iteration with highest lodLimit (checklikelihood_CHR_local.sh)[Top_OM]. 
10.	Combined all files and merged with file (getsnpscall.sh) that matches SNP order in call file to the locus (called POS in the RADSeq files as the unique identifier) [CombQC_NR_snps_callfile,txt] and the file that matches the locus to the SNPID [CombQC_NoReps_SNPIDs.txt] (combine_OM_files.R) [CombQCNR_LM2SC_Map_SNPID].
C.	Combining RADSeq and SNP data
1.	Need to determined the individuals that overlap between the two datasets across the 2 families (Ind_Intersect.R)[RAD_IND_OVERLAP.txt and SNP_IND_OVERLAP.txt]. This is also the file where the replicate individuals in Family 1 was determined [SNP_IND_EXCLUDE.txt]. Order and get names to be the same using vcfR.
2.	Create vcf files for each dataset with only the overlapping individuals using VCFtools (overlap) [Comb_QC_noreps_overlap.recode.vcf and (RAD_2FAMComb_overlap.recode.vcf).
3.	Re-name and order individuals so that they match between the 2 datasets then combine (rename_combine) [RAD_SNP_2FAMcomb.vcf]. 
4.	Create separate family files using VCFtools (by_fam) [RAD_SNP_LM1.recode.vcf and RAD_SNP_LM2.recode.vcf].
5.	Run each family and combined families in LepMap3. Combined families use family 2 during Separate Chromosomes [WithProblems]. 
6.	Problem Markers were removed. Gaps were removed by determining large gaps at end of linkage groups visually, calculating gap length for those linkage groups and removing those markers (analyse_map.R). Markers that were repeatedly miss assigned linkage groups at least 3 times in comparison to: maps for family 1 and family 2, map made from only SNP Chip data, map made from only RAD Seq data and the white spruce map (compare_maps.R). Removed markers with a Cook’s distance of greater than 4N (findoutliers.R). These problematic markers were combined into a single list and removed from a list of only mapped markers. Then the mapped markers without problematic markers were kept from the combined family file [RAD_SNP_2FAMcomb.vcf] using VCFtools (removeSNPs.sh) [RAD_SNP_2FAM_Cook4N.recode.vcf]. 
7.	Map again using LepMap3 on the ARC. 

