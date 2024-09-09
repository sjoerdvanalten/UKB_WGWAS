#We used the following software packages in data analysis: 
#R version 4.0.3 with packages ggplot2, tidyverse and ggExtra; R 4.1.0 with packages qqman, dplyr, plyr, ggplot2, ggrepel, stats, readxl, data.table, 
#fastman; stringr, flexiblas, svMisc, estimatr, testit, doParallel, foreach, bigstatsr, BEDMatrix, lmtest, and AER;
#Plink 1.9 (https://www.cog-genomics.org/plink/1.9); 
#QCTool 2.2 (https://www.well.ox.ac.uk/~gav/qctool_v2/); 
#ldsc 1.01 (https://github.com/bulik/ldsc); 
#BOLT-LMM_v2.4.1 (https://alkesgroup.broadinstitute.org/BOLT-LMM/BOLT-LMM_manual.html)
#python 3.9.5 with packages pandas, numpy, re, os and subprocess.

#The following INPUT files are required to run the code:
#a UKB phenotype basket (e.g., basket2005284.ukb46101_refresh_pheno.csv.gz) with all phenotypes and necessary control variables as listed in preprocess_ukb_v3.py 
#UKB genetic data (imputed, in .bgen forma)
#the list of hapmap SNPs, w_hm3.snplist, can be found here https://zenodo.org/records/7773502)
#The IP weights are called UKBWeightsKFolded.csv in this script. Their construction is described in van Alten, S., Domingue, B. W., Faul, J., Galama, T., & Marees, A. T. (2024). Reweighting UK Biobank 
#+corrects for pervasive selection bias due to volunteering. International Journal of Epidemiology, 53(3), dyae054.
#These weights can be downloaded from the UKB returns catalogue under application ID 55154: https://biobank.ndph.ox.ac.uk/ukb/app.cgi?id=55154

#Run the following files in order:

sbatch Preprocess_ukb.sh #to prepare phenotypes
sbatch GeneticDataUnpack.sh #to unpack genetic data

sbatch GWASOntheWeights.sh #to run the GWAS on the weights (with and without controls)
sbatch GWASOnTheWeightsBolt.sh #to run GWAS on the weights with controls using Bolt-LMM

sbatch WeightedGWASIPW.sh #to run weighted GWAS on the weights
sbatch WeigtedGWASIPWNoControls.sh #to run weighted GWAS on the weights without controls

sbatch GWASIPWAnalyze.sh #to analyze the GWAS on the weights
sbatch GWASIPWNoControlAnalyze.sh #to analyze the GWAS (without controls)
sbatch GWASIPWBoltAnalyze.sh #to analyze GWAS (estimated in Bolt, with controls) 

sbatch GWASOnTheWeightsZoom.sh #to obtain zoomed manhattan plots of top hits

#to obtain summary statistics of our weighted phenotypes:
Rscript WeightedPhenotypes.R

#To run weighted GWAS:
sbatch ResidualizePhenotypes.sh #to residualize these phenotypes from control variables 

#Run GWAS for 10 phenotypes
sbatch WeightedGWASAgeFirstBirth.sh
sbatch WeightedGWASBMI.sh
sbatch WeightedGWASEducation.sh 
sbatch WeightedGWASHeight.sh 
sbatch WeightedGWASSevereObesity.sh 
sbatch WeightedGWASType1Diabetes.sh 
sbatch WeightedGWASHealthRating.sh 
sbatch WeightedGWASBreastCancer.sh
sbatch WeightedGWASPhysicalActivity.sh 
sbatch WeightedGWASDrinksPerWeek.sh 

#Analyze the results to generate the main tables and figures in the paper:
sbatch WeightedGWASAnalyze.sh 

sbatch WeightedGWASSex.sh #Sex is the only GWAS that is not run on any control variable 

#To obtain LD-score regression results (heritability estimates)
sbatch WeightedGWASLDScore.sh 

Rscript GWAS_WGWASCompare.R
Rscript GWASSumSummarize.R 

sbatch GWASSexAnalyze.sh 

sbatch TopHitsZoom.sh 

#Gene-annotation analysis (note: run through fuma.ctglab.nl first
Rscript MAGMAPlots.R 

#Material for supplementary notes:

#Simulation results presented in supplementary note 1:
sbatch SimulateSelection.sh
sbatch SimulateSelectionBinary.sh 
sbatch SimulateSelectionBinary95.sh
Rscript ExtractH2.R 
Rscript ExtractH2Binary.R
Rscript ExctractH2Binary95.R

#compare to top hits to GWAS summary results in the literature
sbatch DetailedLookup.sh

#robustness of the top hits:
sbatch WeightedGWASOneVar.sh
