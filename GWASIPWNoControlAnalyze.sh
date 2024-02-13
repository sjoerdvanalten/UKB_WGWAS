#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 32
#SBATCH -p rome
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/GWASIPWNoControlAnalyze.log 

module load 2021
module load R/4.1.0-foss-2021a


#for c in {1..22}
#do 
#cp ../TEMP/ChromResults/GWAS/IPWNoControls${c}.csv ../TEMP/ChromResults/GWAS/IPWNoControlsChr${c}.tab
#../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/ChromResults/GWAS/IPWNoControlsChr${c}.tab --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field P --out ../TEMP/Clumped/IPWNoControlsClumped${c}
#cp ../TEMP/ChromResults/WGWAS/IPWNoControls${c}.csv ../TEMP/ChromResults/WGWAS/IPWNoControlsChr${c}.tab
#../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/ChromResults/WGWAS/IPWNoControlsChr${c}.tab --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field P --out ../TEMP/Clumped/IPWNoControlsClumped_W${c}
#done 

Rscript WeightedGWASSexAnalyze.R IPWNoControlsChr LassoWeight NA NA NA NA NA NA NA ../TEMP/IPWPhenoNoControls.txt IPWNoControls


module purge 
module load 2021
module load Python/2.7.18-GCCcore-10.3.0-bare

pip install --user numpy
pip install --user bitarray
pip install --user pandas
pip install --user scipy

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/GWAS/IPWNoControls.txt --ignore Z --N-col N_eff --out ../TEMP/IPWGWASNoControls.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/WGWAS/IPWNoControls.txt --ignore Z --N-col N_eff --out ../TEMP/IPWWGWASNoControls.munge --merge-alleles ../INPUT/w_hm3.snplist

../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/IPWGWASNoControls.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/IPWNoControlsResults.h2
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/IPWWGWASNoControls.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/IPWNoControlsWeighted.h2

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/IPWGWASNoControls.munge.sumstats.gz,../TEMP/IPWWGWASNoControls.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/IPWNoControls.rg 