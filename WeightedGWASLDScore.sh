#!/bin/bash
#SBATCH -t 12:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 32
#SBATCH -p thin
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/WeightedGWASLDScore.log 


module purge 
module load 2021
module load Python/2.7.18-GCCcore-10.3.0-bare

pip install --user numpy
pip install --user bitarray
pip install --user pandas
pip install --user scipy

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/TABLES/WeightedMAFFullOutput.txt --p PBootstrap --signed-sumstats ZBootstrap,0 --N 380000 --out ../TEMP/MAFGWAS.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/MAFGWAS.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/MAFGWASResults.h2 

#Heritability estimates of GWAS/WGWAS

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/GWAS/PhysicalActivity.txt --N-col N_eff --signed-sumstats Z,0 --out ../TEMP/PhysicalActivityGWAS.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/WGWAS/PhysicalActivity.txt --N-col N_eff --signed-sumstats Z,0 --out ../TEMP/PhysicalActivityWGWAS.munge --merge-alleles ../INPUT/w_hm3.snplist

../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/PhysicalActivityGWAS.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/PhysicalActivityResults.h2 
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/PhysicalActivityWGWAS.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/PhysicalActivityWeightedResults.h2 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/PhysicalActivityGWAS.munge.sumstats.gz,../TEMP/PhysicalActivityWGWAS.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/PhysicalActivity.rg 


for pheno in {YearsEducation,BMI,Height,DrinksPerWeek,AgeFirstBirth,HealthRating,Type1Diabetes,BreastCancerFemale,SevereObesity}

do 
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/GWAS/${pheno}.txt --N-col N_eff --out ../TEMP/${pheno}GWAS.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/WGWAS/${pheno}.txt --N-col N_eff --out ../TEMP/${pheno}WGWAS.munge --merge-alleles ../INPUT/w_hm3.snplist

#Heritability estimates of GWAS/WGWAS
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/${pheno}GWAS.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/${pheno}Results.h2 
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/${pheno}WGWAS.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/${pheno}WeightedResults.h2 

#Genetic correlations of GWAS/WGWAS
../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/${pheno}GWAS.munge.sumstats.gz,../TEMP/${pheno}WGWAS.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/${pheno}.rg 
done 


	
module load 2021
module load R/4.1.0-foss-2021a

Rscript LDScoreSummarize.R