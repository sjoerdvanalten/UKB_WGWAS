#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 20
#SBATCH -p rome
#SBATCH --output=Test.log 

chr=10
iter=10

#./ldak5.2.linux --bfile ../INPUT/LDAKTestData/hapmap --power 0.25 --her 0.2 --num-phenos 1 --num-causals -1 --make-phenos ../TEMP/TestPheno1

#draw effect sizes from clumped data  


for caus in {1,10,200,2000,-1}
do
module load 2021
module load R/4.1.0-foss-2021a
		
Rscript SimulateSelectionBinary.R $chr $caus $iter 
done