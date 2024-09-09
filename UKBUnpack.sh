#!/bin/bash
#SBATCH -t 1-00:00:00 ## WALL CLOCK TIME
#SBATCH -N 1 -c 32
#SBATCH -p thin
#SBATCH --cpus-per-task=8
#SBATCH --output=UKBUnpack.log 


c
./ukbconv ../GWAS_PIPELINE/INPUT/basket2005284.ukb41217.enc_ukb docs
