#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 20
#SBATCH -p rome
#SBATCH --output=SimulateSelectionBinary.log 

chr=10
iter=30

#./ldak5.2.linux --bfile ../INPUT/LDAKTestData/hapmap --power 0.95 --her 0.2 --num-phenos 1 --num-causals -1 --make-phenos ../TEMP/TestPheno1

#draw effect sizes from clumped data  

../SOFTWARE/plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${chr} --indep-pairwise 500 50 0.5 --threads 10 --out ../TEMP/SIMUL_effect #ensure to remove HLA region chr6:28,477,797-33,448,354 
../SOFTWARE/plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${chr} --extract ../TEMP/SIMUL_effect.prune.in --make-bed --out ../TEMP/SIMUL_effect${chr}


for caus in {1,10,200,2000,-1}
do

./ldak5.2.linux --bfile ../TEMP/SIMUL_effect${chr} --power -0.95 --her 0.2 	--prevalence 0.9 --num-phenos 1 --num-causals ${caus} --make-phenos ../TEMP/TestPhenoBinary95${caus}

#create PGI bsed on simulated effects 
../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef${chr} --score ../TEMP/TestPhenoBinary95${caus}.effects 2 3 6 --out ../TEMP/SimulPGIBinary95${caus}

module load 2021
module load R/4.1.0-foss-2021a
		
Rscript SimulateSelectionBinary95.R $chr $caus $iter 

#run a GWAS 

module load PLINK/2.00a2.3_x86_64


for i in `seq 1 ${iter}` 
do
for j in {1..6}
do 
plink2 --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef10 --keep ../TEMP/SimulSelect1Binary95_${i}${j}_${caus}.txt --linear --pheno ../TEMP/Simul1Binary95_${i}${j}_${caus}.pheno --out ../OUTPUT/Simulations/SimulGWAS_1Binary95${i}${j}_${caus} 
plink2 --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef10 --keep ../TEMP/SimulSelect2aBinary95_${i}${j}_${caus}.txt --linear --pheno ../TEMP/Simul2aBinary95_${i}${j}_${caus}.pheno --out ../OUTPUT/Simulations/SimulGWAS_2aBinary95${i}${j}_${caus} 
plink2 --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef10 --keep ../TEMP/SimulSelect2bBinary95_${i}${j}_${caus}.txt --linear --pheno ../TEMP/Simul2bBinary95_${i}${j}_${caus}.pheno --out ../OUTPUT/Simulations/SimulGWAS_2bBinary95${i}${j}_${caus}
done
done

wait 

for i in `seq 1 ${iter}`
do
for j in {1..6}
do 
Rscript CleanGWASSum.R SimulGWAS_1Binary95${i}${j}_${caus} &
Rscript CleanGWASSum.R SimulGWAS_2aBinary95${i}${j}_${caus} &
Rscript CleanGWASSum.R SimulGWAS_2bBinary95${i}${j}_${caus} &
#../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef10 --clump ../TEMP/SimulGWAS_1${i}${j}_clean.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-snp-field ID --clump-field P --out ../TEMP/SimulClump1${i}${j}
#../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef10 --clump ../TEMP/SimulGWAS_2a${i}${j}_clean.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-snp-field ID --clump-field P --out ../TEMP/SimulClump2a${i}${j}
#../SOFTWARE/./plink --bfile ../TEMP/PLINKFILES/UKBHapMapSNPsDef10 --clump ../TEMP/SimulGWAS_2b${i}${j}_clean.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-snp-field ID --clump-field P --out ../TEMP/SimulClump2b${i}${j}
done
done
wait 

module purge 
module load 2021
module load Python/2.7.18-GCCcore-10.3.0-bare

pip install --user numpy
pip install --user bitarray
pip install --user pandas
pip install --user scipy

#heritability check 

for i in `seq 1 ${iter}`
do
for j in {1..6}
do
../SOFTWARE/ldsc/munge_sumstats.py --sumstats ../TEMP/Simulations/SimulGWAS_1Binary95${i}${j}_${caus}_clean.txt --snp ID --a1 REF --a2 ALT --N-col OBS_CT --p P --out ../TEMP/Simulations/Simul1Binary95${i}${j}_${caus}.munge --merge-alleles ../INPUT/w_hm3.snplist &
../SOFTWARE/ldsc/munge_sumstats.py --sumstats ../TEMP/Simulations/SimulGWAS_2aBinary95${i}${j}_${caus}_clean.txt --snp ID --a1 REF --a2 ALT --N-col OBS_CT --p P --out ../TEMP/Simulations/Simul2aBinary95${i}${j}_${caus}.munge --merge-alleles ../INPUT/w_hm3.snplist &
../SOFTWARE/ldsc/munge_sumstats.py --sumstats ../TEMP/Simulations/SimulGWAS_2bBinary95${i}${j}_${caus}_clean.txt --snp ID --a1 REF --a2 ALT --N-col OBS_CT --p P --out ../TEMP/Simulations/Simul2bBinary95${i}${j}_${caus}.munge --merge-alleles ../INPUT/w_hm3.snplist 

done 
done 

wait

for i in `seq 1 ${iter}`
do
for j in {1..6}
do
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/Simulations/Simul1Binary95${i}${j}_${caus}.munge.sumstats.gz --ref-ld ../INPUT/eur_w_ld_chr/10 --w-ld ../INPUT/eur_w_ld_chr/10 --out ../OUTPUT/Simulations/Simul1Binary95${i}${j}_${caus}.h2 &
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/Simulations/Simul2aBinary95${i}${j}_${caus}.munge.sumstats.gz --ref-ld ../INPUT/eur_w_ld_chr/10 --w-ld ../INPUT/eur_w_ld_chr/10 --out ../OUTPUT/Simulations/Simul2aBinary95${i}${j}_${caus}.h2 &
../SOFTWARE/ldsc/./ldsc.py --h2 ../TEMP/Simulations/Simul2bBinary95${i}${j}_${caus}.munge.sumstats.gz --ref-ld ../INPUT/eur_w_ld_chr/10 --w-ld ../INPUT/eur_w_ld_chr/10 --out ../OUTPUT/Simulations/Simul2bBinary95${i}${j}_${caus}.h2 

done 
done 

wait


rm -r ../TEMP/Simulations/Simul*.munge.sumstats.gz &

#effect sizes check: 
module load 2021
module load R/4.1.0-foss-2021a
	
#Rscript SimulateSelectionEffSizeCheck.R 
done