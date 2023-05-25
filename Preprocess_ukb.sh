#!/bin/bash
#SBATCH -t 12:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 32 ## REQUEST NODES AND CORES (ONLY STAGING NODES CAN REQUEST SINGLE CORES)
#SBATCH --cpus-per-task=32
#SBATCH -p thin ## NODE TYPE: thin, fat, gpu, short, staging
#SBATCH --output=logs/Preprocess_ukb.log 

module load 2021
module load Python/3.9.5-GCCcore-10.3.0-bare 
python preprocess_ukb_v3.py

module load R/4.1.0-foss-2021a
Rscript UKBUnpackControls.R 

