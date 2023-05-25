#!/bin/bash
#SBATCH -t 2-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 128
#SBATCH -p thin
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/WeightedGWASSevereObesity.log 

#Syntax of the arguments is as follows: bim/bed/fam-file, phenotype name, phenotype-file, output file name 
module load 2021
module load R/4.1.0-foss-2021a

TOP5KPath=../../GWAS_SUMMARY/TOP5K

#Rscript UKBSelectGeneticIDs.R
#Note: make sure the name of the Reflist is the same as the name of the residualized phenotypes in Phenotypes.resid.txt 

#for trait in {YearsEducation,BMI,Height,Ashtma,CigsPerDay,SmokingInitiation,Depression,SBP,DBP,WC,WHR,HIP,Drinks}

for c in {1..22}
do
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} SevereObesity ../TEMP/SevereObesity.resid.txt SevereObesityChr${c}.csv ../INPUT/UKBWeightsKFolded.csv
done


#Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} BMI ../TEMP/BMI.resid.txt BMIChr${c}.csv ../INPUT/UKBSelectionWeights.csv

#for trait in {Education,Height,BMI,Ashtma,CigsPerDay,SmokingInitiation,Depression,SBP,DBP,WC,WHR,HIP,Drinks}
#for w in {1..4}
#do 
#Rscript --vanilla WeightedGWASAnalyzeTop5000.R ../TEMP/ResultsLassoYearsEducationVar${w}Chr YearsEducation $TOP5KPath/Top5000Education3.txt "MarkerName" "A1" "A2" "BETA" "SE" $w
#Rscript --vanilla WeightedGWASAnalyzeTop5000.R ../TEMP/ResultsLassoBMIVar${w}Chr BMI $TOP5KPath/Top5000BMI2.txt "SNP" "Tested_Allele" "Other_Allele" "BETA" "SE" $w
#done 