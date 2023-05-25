#!/bin/bash
#SBATCH -t 2-00:00:00 ## WALL CLOCK TIME2
#SBATCH -N 1
#SBATCH -p thin
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/WeightedGWASEducation.log 

#Syntax of the arguments is as follows: bim/bed/fam-file, phenotype name, phenotype-file, output file name 
module load 2021
module load R/4.1.0-foss-2021a

for c in {1..22}
do
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} YearsEducation ../TEMP/YearsEducation.resid.txt YearsEducationChr${c}.tab ../INPUT/UKBWeightsKFolded.csv
done
