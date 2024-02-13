#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 20
#SBATCH -p fat
#SBATCH --output=WeightedGWASOneVar.log 

module load 2021
module load R/4.1.0-foss-2021a
		
Rscript WeightedGWASOneVarT1D.R
Rscript WeightedGWASOneVarBreastCancer.R