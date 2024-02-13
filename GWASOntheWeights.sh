#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1  -c 32
#SBATCH --partition=thin
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/GWASOntheWeights.log 

module load 2021

module load R/4.1.0-foss-2021a

#module load OpenBLAS/0.3.15-GCC-10.3.0 
module load imkl-FFTW/2021.4.0-iimpi-2021b 
#module load ScaLAPACK/2.1.0-gompi-2021a-fb
#module load FlexiBLAS/3.0.4-GCC-10.3.0 
module load NLopt/2.7.0-GCCcore-10.3.0
module load Boost/1.76.0-GCC-10.3.0

Rscript WeightToPheno.R
#without residualizing for covariates

for c in {1..22}
do
#../SOFTWARE/./plink  --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --linear intercept --allow-no-sex --pheno ../TEMP/UKBWeightsPheno.pheno --all-pheno --out ../TEMP/GWASonTheWeights${c}
../SOFTWARE/./plink  --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --linear intercept --allow-no-sex --pheno ../TEMP/UKBWeightsPhenoResid.pheno --all-pheno --out ../TEMP/GWASonTheWeightsControlled${c}
done 



#clump to check for top hits
for c in {1..22}
do
#../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/GWASonTheWeights${c}.P1.assoc.linear  --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-r2 0.1 --clump-field P --out ../TEMP/GWASOnTheWeightsClumped${c}.txt
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/GWASonTheWeightsControlled${c}.P1.assoc.linear  --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-r2 0.1 --clump-field P --out ../TEMP/GWASOnTheWeightsControlledClumped${c}.txt
done 

Rscript GWASOnTheWeightsAnalyze.R ../TEMP/GWASonTheWeights ../TEMP/GWASOnTheWeightsClumped GWASOnTheWeights
Rscript GWASOnTheWeightsAnalyze.R ../TEMP/GWASonTheWeightsControlled ../TEMP/GWASOnTheWeightsControlledClumped GWASOnTheWeightsControlled

module purge 
module load 2021
module load Python/2.7.18-GCCcore-10.3.0-bare

pip install --user numpy
pip install --user bitarray
pip install --user pandas
pip install --user scipy

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../TEMP/GWASOnTheWeightsFull.tab --N-col NMISS --out ../TEMP/GWASOnWeights.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/GWASOnWeights.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/GWASOnWeightsResults.h2 

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../TEMP/GWASOnTheWeightsControlledFull.tab --N-col NMISS --out ../TEMP/GWASOnWeightsControlled.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/GWASOnWeightsControlled.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/GWASOnWeightsControlledResults.h2 

#munge publicly available sumstats:
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/GWAS_EA_excl23andMe.txt  --N 766345 --out ../TEMP/EA3.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/SWB_Full.txt  --N 298000 --out ../TEMP/SWB.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/CigarettesPerDay.txt --a1 ALT --a2 REF --N-col N_EFFECTIVE --out ../TEMP/CigarettesPerDay.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/DrinksPerWeek.txt --a1 REF --a2 ALT --N-col N_EFFECTIVE --out ../TEMP/DrinksPerWeek.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/SmokingInitiation.txt --a1 ALT --a2 REF --N-col N_EFFECTIVE --out ../TEMP/SmokingInitiation.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt --a1 Tested_Allele --a2 Other_Allele --N-col N --out ../TEMP/BMI2.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/GIANT_2015_HIP_COMBINED_EUR.txt --N-col N --out ../TEMP/HIP.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/GIANT_2015_WHR_COMBINED_EUR.txt --N-col N --out ../TEMP/WHR.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/GIANT_2015_WC_COMBINED_EUR.txt --N-col N --out ../TEMP/WC.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/Meta-analysis_Wood_et_al+UKBiobank_2018.txt --a1 Tested_Allele --a2 Other_Allele --N-col N --out ../TEMP/Height.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/daner_PGC_SCZ52_0513a.hq2 --N 150000 --out ../TEMP/Schizophrenia.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/PGC_UKB_depression_genome-wide.txt --signed-sumstats LogOR,0 --N 500000 --out ../TEMP/Depression.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../TEMP/oncoarray_bcac_public_release_oct17_clean.txt.hapmap --ignore Z --N 228951 --a1 a0 --a2 a1 --out ../TEMP/BreastCancer.munge --merge-alleles ../INPUT/w_hm3.snplist
  
awk '!($2="")' ../../svalten/GWAS_SUMMARY/RAW/BP-ICE_EUR_DBP_transformed_15-04-2020.txt  > ../TEMP/BP-ICE_EUR_DBP_transformed_15-04-2020.txt
awk '!($2="")' ../../svalten/GWAS_SUMMARY/RAW/BP-ICE_EUR_SBP_transformed_15-04-2020.txt  > ../TEMP/BP-ICE_EUR_SBP_transformed_15-04-2020.txt
awk '{$2=$9=""; print $0}' ../../svalten/GWAS_SUMMARY/RAW/BP-ICE_EUR_HTN_15-04-2020.txt  > ../TEMP/BP-ICE_EUR_HTN_15-04-2020.txt
awk '!($1="")' ../../svalten/GWAS_SUMMARY/RAW/SHARE-without23andMe.LDSCORE-GC.SE-META.v0   > ../TEMP/SHARE-without23andMe.LDSCORE-GC.SE-META.v0

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../TEMP/BP-ICE_EUR_DBP_transformed_15-04-2020.txt --snp rsID --a1 Allele1 --a2 Allele2 --N-col N --out ../TEMP/DBP.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../TEMP/BP-ICE_EUR_SBP_transformed_15-04-2020.txt --snp rsID --N-col N --out ../TEMP/SBP.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../TEMP/BP-ICE_EUR_HTN_15-04-2020.txt --N-col N --snp rsID --out ../TEMP/HTN.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../TEMP/SHARE-without23andMe.LDSCORE-GC.SE-META.v0 --a1 EFFECT_ALLELE --a2 OTHER_ALLELE --snp RS_ID --N-col N --out ../TEMP/Allergic.munge --merge-alleles ../INPUT/w_hm3.snplist

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/GCST90000048_buildGRCh37.tsv --N 418758 --snp variant_id --out ../TEMP/AFB.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/GCST90012790_buildGRCh38.tsv --N 300639 --snp variant_id --out ../TEMP/ParticipationFoodQuestionnaire.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/GCST90012791_buildGRCh38.tsv --N 215127 --snp variant_id --out ../TEMP/ParticipationPhysicalMonitor.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/GCST90012792_buildGRCh38.tsv --N 294787 --snp variant_id --out ../TEMP/ParticipationMentalHealthQuestionnaire.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../../svalten/GWAS_SUMMARY/RAW/GCST90012793_buildGRCh38.tsv --N 451036 --snp variant_id --out ../TEMP/ParticipationAideMemoire.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../TEMP/GCST90014023_buildGRCh38_cleaned.tsv --N 520580 --p p_value --snp variant_id --out ../TEMP/Type1Diabetes.munge --merge-alleles ../INPUT/w_hm3.snplist


../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/EA3.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandEA3.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/CigarettesPerDay.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandCigarettesPerDay.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/DrinksPerWeek.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandDrinksPerWeek.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/SmokingInitiation.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandSmokingInitiation.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/BMI2.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandBMI2.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/HIP.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandHIP.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/WHR.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandWHR.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/WC.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandWC.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/Height.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandHeight.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/Schizophrenia.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandSchizophrenia.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/Depression.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandDepression.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/DBP.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandDBP.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/SBP.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandSBP.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/HTN.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandHTN.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/Allergic.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandAllergic.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/SWB.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandSWB.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/AFB.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandAFB.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/ParticipationFoodQuestionnaire.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandParticipationFoodQuestionnaire.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/ParticipationPhysicalMonitor.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandParticipationPhysicalMonitor.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/ParticipationMentalHealthQuestionnaire.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandParticipationMentalHealthQuestionnaire.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/ParticipationAideMemoire.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandParticipationAideMemoire.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/Type1Diabetes.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandType1Diabetes.rg 

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/BreastCancer.munge.sumstats.gz,../TEMP/GWASOnWeights.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/GWASWeightsandBreastCancer.rg 

module purge 
module load 2021
module load R/4.1.0-foss-2021a

Rscript GeneticCorrelationsGWASOnTheWeights.R 


#Rscript MAFBootstrapTopHitsPlot.R	

#MAF values in LDScore pipeline