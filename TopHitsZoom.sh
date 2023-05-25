#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1  -c 16
#SBATCH --partition=fat
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/TopHitsZoom.log 

module load 2021

module load R/4.1.0-foss-2021a

#module load OpenBLAS/0.3.15-GCC-10.3.0 
module load imkl-FFTW/2021.4.0-iimpi-2021b 
#module load ScaLAPACK/2.1.0-gompi-2021a-fb
#module load FlexiBLAS/3.0.4-GCC-10.3.0 
module load NLopt/2.7.0-GCCcore-10.3.0
module load Boost/1.76.0-GCC-10.3.0

#zoom in on newly found loci in WGWAS. (so far, we have found these for breast cancer and T1D).
Rscript UKBTopHitsRanges.R

for c in {1,3,4,5,6,7,11,12,14,15,17,18} #only include chromosomes for which tophits were found 
do

../SOFTWARE/qctool/build/release/apps/./qctool_v2.2.0 -g /projects/0/koelling/data/UKB/UKB_v3/imp/ukb_imp_chr${c}_v3.bgen -s /projects/0/koelling/data/UKB/UKB_v3/imp/sample/ukb_imp_chr1_22_v3.sample -incl-ranges ../TEMP/TopHitsBPRangeHausman.txt -ofiletype binary_ped -og ../TEMP/PLINKFILES/UKBAroundTophits${c}

sort ../TEMP/PLINKFILES/UKBAroundTophits${c}.bim | awk '{print$2}' |uniq -d > ../TEMP/UKBAroundTophits${c}.remove
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBAroundTophits${c} --exclude ../TEMP/UKBAroundTophits${c}.remove --keep ../INPUT/GeneIDs.txt --geno 0.02 --maf 0.01 --hwe 10e-6 --make-bed --out ../TEMP/PLINKFILES/UKBAroundTophitsQC${c}

#calculate r2 around each tophit 
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBAroundTophitsQC${c} --r2 --ld-snp-list ../TEMP/TopHitsSNPRangeHausman.txt --ld-window 5000 --ld-window-kb 500 --ld-window-r2 0.1 --out ../TEMP/PLINKFILES/UKBAroundTophitsQC${c} 

awk 'NR!=1{print$6}' ../TEMP/PLINKFILES/UKBAroundTophitsQC${c}.ld > ../TEMP/TopHits${c}LDSNPs.txt
#Redo the weighted MAF analyses on the regions surrounding top hits only (+/- 500 kb), in order to get a good picture of these regions, using not only HapMap SNPs, but all SNPs present in UKB
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBAroundTophitsQC${c} --snps-only --extract ../TEMP/TopHits${c}LDSNPs.txt --make-bed --out ../TEMP/PLINKFILES/UKBAroundTophitsQC${c}LD0.1 

done 

Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC1LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr1.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC3LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr3.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC5LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr5.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC6LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr6.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC7LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr7.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC11LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr11.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC12LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr12.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC14LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr14.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC15LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr15.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC17LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr17.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC18LD0.1 Type1Diabetes ../TEMP/Type1Diabetes.resid.txt HausmanAroundTophitsChr18.csv ../INPUT/UKBWeightsKFolded.csv

Rscript UKBTopHitsPlot.R ../OUTPUT/TABLES/DifferentSNPsWGWASSuggestiveType1Diabetes.csv

Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC3LD0.1 BreastCancer ../TEMP/BreastCancer.resid.txt HausmanAroundTophitsChr3.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC4LD0.1 BreastCancer ../TEMP/BreastCancer.resid.txt HausmanAroundTophitsChr4.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC5LD0.1 BreastCancer ../TEMP/BreastCancer.resid.txt HausmanAroundTophitsChr5.csv ../INPUT/UKBWeightsKFolded.csv
Rscript --vanilla WeightedGWAS.R ../TEMP/PLINKFILES/UKBAroundTophitsQC11LD0.1 BreastCancer ../TEMP/BreastCancer.resid.txt HausmanAroundTophitsChr11.csv ../INPUT/UKBWeightsKFolded.csv

Rscript UKBTopHitsPlot.R ../OUTPUT/TABLES/DifferentSNPsWGWASSuggestiveBreastCancer.csv
