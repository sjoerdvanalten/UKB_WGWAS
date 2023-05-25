#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 32
#SBATCH -p thin
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/WeightedGWASAnalyze.log 


module load 2021
module load R/4.1.0-foss-2021a

#Calculate MAFs in data:

for c in {1..22}
do
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --freq --out ../TEMP/MAF${c}
done

#Clump original GWAS summary statistics (when they exist)
for c in {1..22} 
do
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../../svalten/GWAS_SUMMARY/SNP_gwas_mc_merge_nogc.tbl.uniq.hapmap --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field p  --clump-r2 0.1 --out ../TEMP/GWASSum/BMI1Chr${c}
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../../svalten/GWAS_SUMMARY/EduYears_Main.txt.hapmap --clump-p1 1 --clump-p2 1 --clump-snp-field MarkerName --clump-field Pval  --clump-r2 0.1 --out ../TEMP/GWASSum/EA2Chr${c}
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../../svalten/GWAS_SUMMARY/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.hapmap --clump-p1 1 --clump-p2 1  --clump-r2 0.1 --clump-snp-field MarkerName --clump-field p --out ../TEMP/GWASSum/Height1Chr${c}
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../../svalten/GWAS_SUMMARY/DrinksPerWeek.WithoutUKB.txt.hapmap  --clump-p1 1 --clump-p2 1  --clump-r2 0.1 --clump-snp-field RSID --clump-field PVALUE --out ../TEMP/GWASSum/DrinksPerWeekChr${c}
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/oncoarray_bcac_public_release_oct17_clean.txt.hapmap  --clump-p1 1 --clump-p2 1  --clump-r2 0.1 --clump-snp-field SNP --clump-field P --out ../TEMP/GWASSum/BreastCancerChr${c}
done


#Clump GWAS/WGWAS summary statitistics 
for c in {1..22}
do 
for pheno in {YearsEducation,BMI,Height,AgeFirstBirth,SevereObesity,BreastCancer,PhysicalActivity,HealthRating,Type1Diabetes,DrinksPerWeek}
do
echo ${c} ${pheno}
Rscript CSVToTab.R ../TEMP/ChromResults/GWAS/${pheno}Chr${c}.csv ../TEMP/ChromResults/GWAS/${pheno}Chr${c}.tab 
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/ChromResults/GWAS/${pheno}Chr${c}.tab --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field P --out ../TEMP/Clumped/${pheno}Clumped${c}
Rscript CSVToTab.R ../TEMP/ChromResults/WGWAS/${pheno}Chr${c}.csv ../TEMP/ChromResults/WGWAS/${pheno}Chr${c}.tab 
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/ChromResults/WGWAS/${pheno}Chr${c}.tab --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field P --out ../TEMP/Clumped/${pheno}Clumped_W${c}
done 
done

Rscript WeightedGWASAnalyze.R YearsEducationChr YearsEducation ../../svalten/GWAS_SUMMARY/EduYears_Main.txt.hapmap MarkerName A1 A2 BETA SE ../TEMP/GWASSum/EA2Chr ../TEMP/YearsEducation.resid.txt YearsEducation
Rscript WeightedGWASAnalyze.R BMIChr BMI ../../svalten/GWAS_SUMMARY/SNP_gwas_mc_merge_nogc.tbl.uniq.hapmap SNP A1 A2 BETA SE ../TEMP/GWASSum/BMI1Chr ../TEMP/BMI.resid.txt BMI
Rscript WeightedGWASAnalyze.R HeightChr Height ../../svalten/GWAS_SUMMARY/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.hapmap MarkerName Allele1 Allele2 BETA SE ../TEMP/GWASSum/Height1Chr ../TEMP/Height.resid.txt Height
Rscript WeightedGWASAnalyze.R DrinksPerWeekChr DrinksPerWeek ../../svalten/GWAS_SUMMARY/DrinksPerWeek.WithoutUKB.txt.hapmap RSID REF ALT BETA SE ../TEMP/GWASSum/DrinksPerWeekChr ../TEMP/DrinksPerWeek.resid.txt DrinksPerWeek
Rscript WeightedGWASAnalyze.R BreastCancerFemaleChr BreastCancer ../TEMP/oncoarray_bcac_public_release_oct17_clean.txt.hapmap SNP a0 a1 BETA SE ../TEMP/GWASSum/BreastCancerChr ../TEMP/BreastCancer.Female.resid.txt BreastCancer

Rscript WeightedGWASAnalyze.R SevereObesityChr SevereObesity ../../svalten/GWAS_SUMMARY/SNP_gwas_mc_merge_nogc.tbl.uniq.hapmap SNP A1 A2 BETA SE ../TEMP/GWASSum/BMI1Chr ../TEMP/SevereObesity.resid.txt SevereObesity

Rscript WeightedGWASAnalyze.R AgeFirstBirthChr AgeFirstBirth NA NA NA NA NA NA NA ../TEMP/AgeFirstBirth.resid.txt AgeFirstBirth
Rscript WeightedGWASAnalyze.R PhysicalActivityChr PhysicalActivity NA NA NA NA NA NA NA ../TEMP/PhysicalActivity.resid.txt PhysicalActivity
Rscript WeightedGWASAnalyze.R HealthRatingChr HealthRating NA NA NA NA NA NA NA ../TEMP/HealthRating.resid.txt HealthRating
Rscript WeightedGWASAnalyze.R Type1DiabetesChr Type1Diabetes NA NA NA NA NA NA NA ../TEMP/Type1Diabetes.resid.txt Type1Diabetes



