#!/bin/bash
#SBATCH -t 4-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 128 ## REQUEST NODES AND CORES (ONLY STAGING NODES CAN REQUEST SINGLE CORES)
#SBATCH -p fat_rome ## NODE TYPE: thin, fat, gpu, short, staging
#SBATCH --cpus-per-task=32
#SBATCH --output=logs/GeneticDataUnpack.log 

awk 'NR!=1{print$1}' ../INPUT/w_hm3.snplist > ../TEMP/hm3_rsid.txt
awk 'NR!=1{print $1,$2}' ../INPUT/w_hm3.snplist  > ../TEMP/hm3reflist.txt
#c=$SLURM_ARRAY_TASK_ID

#Syntax of the arguments is as follows: bim/bed/fam-file, phenotype name, phenotype-file, output file name 
module load 2021	
#module load PLINK/2.00a2.3_x86_64	
module load foss/2021a
module load Boost/1.76.0-GCC-10.3.0
module load Eigen/3.3.9-GCCcore-10.3.0
module load SQLite/3.35.4-GCCcore-10.3.0
module load Python/3.9.5-GCCcore-10.3.0

for c in `seq 1 22`
 do
 #../SOFTWARE/qctool/build/release/apps/./qctool_v2.2.0 -g /projects/0/koelling/data/UKB/UKB_v3/imp/ukb_imp_chr${c}_v3.bgen -s /projects/0/koelling/data/UKB/UKB_v3/imp/sample/ukb_imp_chr1_22_v3.sample -incl-rsids ../TEMP/hm3_rsid.txt -ofiletype binary_ped -og ../TEMP/PLINKFILES/UKBHapMapSNPs${c}  
 ../SOFTWARE/qctool/build/release/apps/./qctool_v2.2.0 -g /projects/0/koelling/data/UKB/UKB_v3/imp/ukb_imp_chr${c}_v3.bgen -s /projects/0/koelling/data/UKB/UKB_v3/imp/sample/ukb_imp_chr1_22_v3.sample -ofiletype binary_ped -og ../TEMP/PLINKFILES/UKBAllSNPs${c}  
 sort ../TEMP/PLINKFILES/UKBAllSNPs${c}.bim | awk '{print$2}' |uniq -d > ../TEMP/UKBAllSNPs${c}.remove
 ../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBAllSNPs${c} --keep ../INPUT/GeneIDs.txt --geno 0.02 --maf 0.01 --hwe 10e-6 --make-bed --out ../TEMP/PLINKFILES/UKBAllSNPsDef${c}
 rm -r ../TEMP/PLINKFILES/UKBAllSNPs1.bed ../TEMP/PLINKFILES/UKBAllSNPs1.bim ../TEMP/PLINKFILES/UKBAllSNPs1.fam
 gzip ../TEMP/PLINKFILES/UKBAllSNPsDef${c}.bed
 ../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPs${c} --exclude ../TEMP/UKBAllSNPs${c}.remove --keep ../INPUT/GeneIDs.txt --geno 0.02 --maf 0.01 --hwe 10e-6 --reference-allele ../TEMP/hm3reflist.txt --make-bed --out ../TEMP/PLINKFILES/UKBHapMapSNPsDef${c}&	
done
wait

#filter multiallelic SNPs and exclude them in PLINK: (to do?)