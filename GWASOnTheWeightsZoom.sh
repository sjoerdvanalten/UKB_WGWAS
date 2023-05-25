#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1  -c 32
#SBATCH --partition=thin
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/GWASOntheWeightsZoom.log 

module load 2021

module load R/4.1.0-foss-2021a

#module load OpenBLAS/0.3.15-GCC-10.3.0 
module load imkl-FFTW/2021.4.0-iimpi-2021b 
#module load ScaLAPACK/2.1.0-gompi-2021a-fb
#module load FlexiBLAS/3.0.4-GCC-10.3.0 
module load NLopt/2.7.0-GCCcore-10.3.0
module load Boost/1.76.0-GCC-10.3.0

for c in {1,2,4,6,13} #only include chromosomes for which tophits were found 
do

../SOFTWARE/qctool/build/release/apps/./qctool_v2.2.0 -g /projects/0/koelling/data/UKB/UKB_v3/imp/ukb_imp_chr${c}_v3.bgen -s /projects/0/koelling/data/UKB/UKB_v3/imp/sample/ukb_imp_chr1_22_v3.sample -incl-ranges ../TEMP/TopHitsBPRangeGWASOnTheWeights.txt -ofiletype binary_ped -og ../TEMP/PLINKFILES/UKBGWASAroundTophits${c}

sort ../TEMP/PLINKFILES/UKBGWASAroundTophits${c}.bim | awk '{print$2}' |uniq -d > ../TEMP/UKBGWASAroundTophits${c}.remove
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBGWASAroundTophits${c} --exclude ../TEMP/UKBGWASAroundTophits${c}.remove --keep ../INPUT/GeneIDs.txt --geno 0.02 --maf 0.01 --hwe 10e-6 --make-bed --out ../TEMP/PLINKFILES/UKBGWASAroundTophitsQC${c}

#calculate r2 around each tophit 
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBGWASAroundTophitsQC${c} --r2 --ld-snp-list ../TEMP/TopHitsSNPRangeGWASOnTheWeights.txt --ld-window 5000 --ld-window-kb 500 --ld-window-r2 0.1 --out ../TEMP/PLINKFILES/UKBGWASAroundTophitsQC${c} 

awk 'NR!=1{print$6}' ../TEMP/PLINKFILES/UKBGWASAroundTophitsQC${c}.ld > ../TEMP/TopHits${c}LDSNPs.txt
#Redo the weighted MAF analyses on the regions surrounding top hits only (+/- 500 kb), in order to get a good picture of these regions, using not only HapMap SNPs, but all SNPs present in UKB
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBGWASAroundTophitsQC${c} --extract ../TEMP/TopHits${c}LDSNPs.txt --make-bed --out ../TEMP/PLINKFILES/UKBGWASAroundTophitsQC${c}LD0.1 

../SOFTWARE/./plink  --bfile ../TEMP/PLINKFILES/UKBGWASAroundTophitsQC${c}LD0.1 --linear --allow-no-sex --pheno ../TEMP/UKBWeightsPheno.pheno --all-pheno --out ../TEMP/GWASonTheWeightsAroundTopHits${c}
Rscript --vanilla MAFsBootstrapPValues.R ../TEMP/PLINKFILES/UKBAroundTophitsQC${c}LD0.1 TopHits${c} 16
done 

Rscript GWASOnTheWeightsTopHitsPlot.R
