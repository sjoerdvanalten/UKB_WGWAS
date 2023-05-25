#!/bin/bash
#SBATCH -t 12:00:00 ## WALL CLOCK TIME
#SBATCH -n 1
#SBATCH -N 1 -c 20
#SBATCH -p thin 
#SBATCH --output=logs/ResidualizePhenotypes.log 

INPUTDIR=../INPUT/PHENOTYPES/

module load 2021

module load R/4.1.0-foss-2021a

Rscript UKBUnpackControls.R 

Rscript --vanilla ResidualizePhenotypes.R ${INPUTDIR}/ageFirstBirth_pheno.PREPARED.txt AgeFirstBirth
Rscript --vanilla ResidualizePhenotypes.R ${INPUTDIR}/bmi_pheno.PREPARED.txt BMI 
Rscript --vanilla ResidualizePhenotypes.R ${INPUTDIR}/cancerBreast_pheno.PREPARED.txt BreastCancer
Rscript --vanilla ResidualizePhenotypes.R ${INPUTDIR}/t1d_pheno.PREPARED.txt Type1Diabetes
Rscript --vanilla ResidualizePhenotypes.R ${INPUTDIR}/height_pheno.PREPARED.txt Height
Rscript --vanilla ResidualizePhenotypes.R ${INPUTDIR}/dpw_pheno.PREPARED.txt DrinksPerWeek
Rscript --vanilla ResidualizePhenotypes.R ${INPUTDIR}/obesitySevere_pheno.PREPARED.txt SevereObesity     	
Rscript --vanilla ResidualizePhenotypes.R ${INPUTDIR}/educYears_pheno.PREPARED.txt YearsEducation            
Rscript --vanilla ResidualizePhenotypes.R ${INPUTDIR}/healthRating_pheno.PREPARED.txt HealthRating        
Rscript --vanilla ResidualizePhenotypes.R ${INPUTDIR}/actModVig_pheno.PREPARED.txt PhysicalActivity