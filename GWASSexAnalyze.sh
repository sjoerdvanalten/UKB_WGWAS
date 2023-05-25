#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 32
#SBATCH -p thin
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/GWASSexAnalyze.log 

module load 2021
module load R/4.1.0-foss-2021a


for c in {1..22}
do 
cp ../TEMP/ChromResults/GWAS/Sex${c}.csv ../TEMP/ChromResults/GWAS/SexChr${c}.tab
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/ChromResults/GWAS/SexChr${c}.tab --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field P --out ../TEMP/Clumped/SexClumped${c}
cp ../TEMP/ChromResults/WGWAS/Sex${c}.csv ../TEMP/ChromResults/WGWAS/SexChr${c}.tab
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c} --clump ../TEMP/ChromResults/WGWAS/SexChr${c}.tab --clump-p1 1 --clump-p2 1 --clump-snp-field SNP --clump-field P --out ../TEMP/Clumped/SexClumped_W${c}
done 

Rscript WeightedGWASSexAnalyze.R SexChr Sex NA NA NA NA NA NA NA ../TEMP/sex.txt Sex


module purge 
module load 2021
module load Python/2.7.18-GCCcore-10.3.0-bare

pip install --user numpy
pip install --user bitarray
pip install --user pandas
pip install --user scipy

../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/GWAS/Sex.txt --ignore Z --N-col N_eff --out ../TEMP/SexGWAS.munge --merge-alleles ../INPUT/w_hm3.snplist
../SOFTWARE/ldsc/./munge_sumstats.py --sumstats ../OUTPUT/WGWAS/Sex.txt --ignore Z --N-col N_eff --out ../TEMP/SexWGWAS.munge --merge-alleles ../INPUT/w_hm3.snplist

../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/SexGWAS.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/SexResults.h2
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/SexWGWAS.munge.sumstats.gz --ref-ld-chr ../INPUT/eur_w_ld_chr/ --w-ld-chr ../INPUT/eur_w_ld_chr/ --out ../OUTPUT/SexWeighted.h2

../SOFTWARE/ldsc/./ldsc.py \
--rg ../TEMP/SexGWAS.munge.sumstats.gz,../TEMP/SexWGWAS.munge.sumstats.gz  \
--ref-ld-chr ../INPUT/eur_w_ld_chr/ \
--w-ld-chr ../INPUT/eur_w_ld_chr/ \
--out ../OUTPUT/Sex.rg 