#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 32
#SBATCH -p rome
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/GWASIPWBoltAnalyze.log 

module load 2021
module load R/4.1.0-foss-2021a


#for c in {1..22}
#do 
#cp ../TEMP/ChromResults/GWAS/IPWNoControls${c}.csv ../TEMP/ChromResults/GWAS/IPWNoControlsChr${c}.tab
#../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/ChromResults/GWAS/IPWNoControlsChr${c}.tab --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field P --out ../TEMP/Clumped/IPWNoControlsClumped${c}
#cp ../TEMP/ChromResults/WGWAS/IPWNoControls${c}.csv ../TEMP/ChromResults/WGWAS/IPWNoControlsChr${c}.tab
#../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/ChromResults/WGWAS/IPWNoControlsChr${c}.tab --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field P --out ../TEMP/Clumped/IPWNoControlsClumped_W${c}
#done 

#Rscript WeightedGWASSexAnalyze.R IPWNoControlsChr LassoWeight NA NA NA NA NA NA NA ../TEMP/IPWPhenoNoControls.txt IPWNoControls


module purge 
module load 2021
module load Python/2.7.18-GCCcore-10.3.0-bare

pip install --user numpy
pip install --user bitarray
pip install --user pandas
pip install --user scipy

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/BOLT.PLINK.Weights.Bgen --ignore Z --N 380000 --a1 ALLELE1 --a2 ALLELE0 --p P_BOLT_LMM_INF --out ../TEMP/IPWGWASBolt.munge --merge-alleles ../INPUT/w_hm3.snplist

../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/IPWGWASBolt.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/IPWGWASBoltResults.h2

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/IPWGWAS.munge.sumstats.gz,../TEMP/IPWGWASBolt.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/Bolt.rg 