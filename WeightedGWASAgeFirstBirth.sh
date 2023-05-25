#!/bin/bash
#SBATCH -t 2-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 128
#SBATCH -p thin
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/WeightedGWASAgeFirstBirth.log 

#Syntax of the arguments is as follows: bim/bed/fam-file, phenotype name, phenotype-file, output file name 
module load 2021
module load R/4.1.0-foss-2021a

for c in {1..22}
do
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} AgeFirstBirth ../TEMP/AgeFirstBirth.resid.txt AgeFirstBirthChr${c}.csv ../INPUT/UKBWeightsKFolded.csv
done

