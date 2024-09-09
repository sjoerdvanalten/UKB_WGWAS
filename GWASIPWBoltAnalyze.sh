#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 32
#SBATCH -p rome
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/GWASIPWBoltAnalyze.log 

module purge 
module load 2021
module load Python/2.7.18-GCCcore-10.3.0-bare

pip install --user numpy
pip install --user bitarray
pip install --user pandas
pip install --user scipy

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/GWAS/IPWBolt.txt --ignore Z --N 380000 --a1 ALLELE1 --a2 ALLELE0 --p P_BOLT_LMM_INF --out ../TEMP/IPWGWASBolt.munge --merge-alleles ../INPUT/w_hm3.snplist

../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/IPWGWASBolt.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/IPWGWASBoltResults.h2

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/IPWGWAS.munge.sumstats.gz,../TEMP/IPWGWASBolt.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/Bolt.rg 