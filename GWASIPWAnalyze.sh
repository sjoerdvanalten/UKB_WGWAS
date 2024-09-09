#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 32
#SBATCH -p rome
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/GWASIPWAnalyze.log 

module load 2021
module load R/4.1.0-foss-2021a


for c in {1..22}
do 
cp ../TEMP/ChromResults/GWAS/IPW${c}.csv ../TEMP/ChromResults/GWAS/IPWChr${c}.tab
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/ChromResults/GWAS/IPWChr${c}.tab --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field P --out ../TEMP/Clumped/IPWClumped${c}
cp ../TEMP/ChromResults/WGWAS/IPW${c}.csv ../TEMP/ChromResults/WGWAS/IPWChr${c}.tab
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/ChromResults/WGWAS/IPWChr${c}.tab --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field P --out ../TEMP/Clumped/IPWClumped_W${c}
done 

Rscript WeightedGWASSexAnalyze.R IPWChr LassoWeight NA NA NA NA NA NA NA ../TEMP/IPWPheno.txt IPW


module purge 
module load 2021
module load Python/2.7.18-GCCcore-10.3.0-bare

pip install --user numpy
pip install --user bitarray
pip install --user pandas
pip install --user scipy

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/GWAS/IPW.txt --ignore Z --N-col N_eff --out ../TEMP/IPWGWAS.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/WGWAS/IPW.txt --ignore Z --N-col N_eff --out ../TEMP/IPWWGWAS.munge --merge-alleles ../INPUT/w_hm3.snplist

../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/IPWGWAS.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/IPWResults.h2
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/IPWWGWAS.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/IPWWeighted.h2

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/IPWGWAS.munge.sumstats.gz,../TEMP/IPWWGWAS.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/IPW.rg 
