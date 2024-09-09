#!/bin/bash
#SBATCH -t 5-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 
#SBATCH --partition=fat_rome
#SBATCH --output=logs/GWASOntheWeightsBolt.log 
#SBATCH --ntasks=1

module load 2021

module load R/4.1.0-foss-2021a

Rscript WeightToPheno.R

module load foss/2021a
module load Boost/1.76.0-GCC-10.3.0
module load Eigen/3.3.9-GCCcore-10.3.0
module load SQLite/3.35.4-GCCcore-10.3.0
module load Python/3.9.5-GCCcore-10.3.0

sort ../INPUT/UKBv3_maf01_info60_rsq30.bim | awk '{print$2}' |uniq -d > ../TEMP/BoltSNPs.txt

#create bolt sample file file 
awk '{print $1, $2}' ../TEMP/PLINKFILES/UKBHapMapSNPsDef22.fam  > ../TEMP/samples.keep
awk '{print $2, $3, $4, $5, $6}' ../TEMP/PLINKFILES/UKBHapMapSNPsDef22.fam  > ../TEMP/PLINKFILES/BoltSNPs22.fam
cp ../TEMP/PLINKFILES/UKBHapMapSNPsDef22.bim tmp.bim 
cp ../TEMP/PLINKFILES/UKBHapMapSNPsDef22.bed tmp.bed 
cp ../TEMP/PLINKFILES/UKBHapMapSNPsDef22.fam tmp.fam 

../SOFTWARE/./plink --bfile tmp --make-bed --out ../TEMP/PLINKFILES/BoltSNPsDef22

awk '{print $2, $2, $3, $4, $5, $6}' ../TEMP/PLINKFILES/BoltSNPsDef22.fam  > ../TEMP/PLINKFILES/bolt.fam
awk 'BEGIN {print "0 0"} {print $0}' tmp > tmp2
awk 'BEGIN {print "ID_1 ID_2"} {print $0}' tmp2 > ../TEMP/samplesbolt.sample
awk '{print $1, $2}' ../TEMP/samples.keep  > tmp2
awk '{print $1}' ../TEMP/PLINKFILES/bolt.fam  > ../TEMP/SampleList.txt

#Extract pruned SNPs 
for c in `seq 1 22`
do
 ../SOFTWARE/qctool/build/release/apps/./qctool_v2.2.0 -g /projects/0/koelling/data/UKB/UKB_v3/imp/ukb_imp_chr${c}_v3.bgen -s /projects/0/koelling/data/UKB/UKB_v3/imp/sample/ukb_imp_chr1_22_v3.sample -incl-rsids ../TEMP/BoltSNPs.txt -incl-samples ../TEMP/SampleList.txt -ofiletype binary_ped -og ../TEMP/PLINKFILES/BoltSNPs${c} &
done 
wait 

# Update FID to be the same as IID
for c in `seq 1 22`
do
  # Create a temporary file with FID=IID
  awk '{print $1, $2, $2, $2}' ../TEMP/PLINKFILES/BoltSNPs${c}.fam > ../TEMP/PLINKFILES/update_ids.txt
  
  # Update FID in the PLINK binary files
  ../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/BoltSNPs${c} --update-ids ../TEMP/PLINKFILES/update_ids.txt --make-bed --out ../TEMP/PLINKFILES/BoltSNPs${c}_updated
  
  # Replace the original files with the updated ones
  mv ../TEMP/PLINKFILES/BoltSNPs${c}_updated.bed ../TEMP/PLINKFILES/BoltSNPs${c}.bed
  mv ../TEMP/PLINKFILES/BoltSNPs${c}_updated.bim ../TEMP/PLINKFILES/BoltSNPs${c}.bim
  mv ../TEMP/PLINKFILES/BoltSNPs${c}_updated.fam ../TEMP/PLINKFILES/BoltSNPs${c}.fam
done

module purge 

module load 2023
module load Boost/1.82.0-GCC-12.3.0
module load NLopt/2.7.1-GCCcore-12.3.0
module load OpenMPI/4.1.5-GCC-12.3.0
module load intel/2023a 

echo "starting BOLT"

# Get the number of CPU cores
AVAILABLE_CORES=$(nproc)

# Subtract 1 from the total
cores=$((AVAILABLE_CORES - 1))

head ../TEMP/samplesbolt.sample
## Note: SPECIFY ONLY 1 FAM FILE ##
../SOFTWARE/BOLT-LMM_v2.4.1/bolt --lmm \
   --LDscoresFile=BOLT-LMM_v2.4.1/tables/LDSCORE.1000G_EUR.tab.gz \
  --bgenMinMAF=0.01 \
  --bgenMinINFO=0.7 \
  --fam=../TEMP/PLINKFILES/BoltSNPs1.fam \
  --bed=../TEMP/PLINKFILES/BoltSNPs{1..22}.bed \
  --bim=../TEMP/PLINKFILES/BoltSNPs{1..22}.bim \
  --phenoFile=../TEMP/UKBWeightsPhenoResid.pheno_bolt \
  --phenoCol=LassoWeightResid \
  --statsFile=../OUTPUT/BOLT.PLINK.Weights \
  --bgenFile=/projects/0/koelling/data/UKB/UKB_v3/imp/ukb_imp_chr{1..22}_v3.bgen \
  --sampleFile=/projects/0/koelling/data/UKB/UKB_v3/imp/sample/ukb_imp_chr1_22_v3.sample  \
  --statsFileBgenSnps=../OUTPUT/BOLT.PLINK.Weights.Bgen \
  --numThreads $cores 
wait

cut -f1 ../INPUT/w_hm3.snplist > tmp.txt
grep -wFf tmp.txt ../OUTPUT/BOLT.PLINK.Weights.Bgen | awk '{$1=$1; print}' OFS=" " > ../OUTPUT/GWAS/IPWBolt.txt
