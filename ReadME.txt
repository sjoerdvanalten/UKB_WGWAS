We used the following software packages in data analysis: 
R version 4.0.3 with packages ggplot2, tidyverse and ggExtra; R 4.1.0 with packages qqman, dplyr, plyr, ggplot2, ggrepel, stats, readxl, data.table, 
fastman; stringr, flexiblas, svMisc, estimatr, testit, doParallel, foreach, bigstatsr, BEDMatrix, lmtest, and AER;
Plink 1.9 (https://www.cog-genomics.org/plink/1.9); 
QCTool 2.2 (https://www.well.ox.ac.uk/~gav/qctool_v2/); 
ldsc 1.01 (https://github.com/bulik/ldsc); 
python 3.9.5 with packages pandas, numpy, re, os and subprocess.

#Run the following files in order:

sbatch Preprocess_ukb.sh to #prepare phenotypes
sbatch GeneticDataUnpack.sh #to unpack genetic data

sbatch GWASOntheWeights.sh #to run the GWAS on the weights (with and without controls)
sbatch WeightedGWASIPW.sh #to run weighted GWAS on the weights
sbatch WeigtedGWASIPWNoControls.sh #to run weighted GWAS on the weights without controls

sbatch GWASIPWAnalyze.sh
sbatch GWASIPWNoControlAnalyze.sh 
sbatch GWASOnTheWeightsZoom.sh to obtain zoomed manhattan plots of top hits

#to obtain summary statistics of our weighted phenotypes:
Rscript WeightedPhenotypes.R

#To run weighted GWAS:
sbatch ResidualizePhenotypes.sh to residualize these phenotypes from control variables 

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
