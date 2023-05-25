#Run the following files in order:

Preprocess_ukb.sh to prepare phenotypes
GeneticDataUnpack.sh to unpack genetic data

GWASOntheWeights.sh to run the GWAS on the weights
GWASOnTheWeightsZoom.sh to obtain zoomed manhattan plots of top hits

#to obtain summary statistics of our weighted phenotypes:
WeightedPhenotypes.R

#To run weighted GWAS:
ResidualizePhenotypes.sh to residualize these phenotypes from control variables 

#Run GWAS for 10 phenotypes
WeightedGWASAgeFirstBirth.sh
WeightedGWASBMI.sh
WeightedGWASEducation.sh 
WeightedGWASHeight.sh 
WeightedGWASSevereObesity.sh 
WeightedGWASType1Diabetes.sh 
WeightedGWASHealthRating.sh 
WeightedGWASBreastCancer.sh
WeightedGWASPhysicalActivity.sh 
WeightedGWASDrinksPerWeek.sh 

WeightedGWASAnalyze.sh


WeightedGWASSex.sh #Sex is the only GWAS that is not run on any control variable 



#To obtain LD-score regression results (heritability estimates)
WeightedGWASLDScore.sh 

Rscript GWAS_WGWASCompare.R
Rscript GWASSumSummarize.R 

GWASSexAnalyze.sh 

TopHitsZoom.sh 

#Gene-annotation analysis (note: run through fuma.ctglab.nl first
MAGMAPlots.R 

#Simulation results presented in supplementary note 1:

SelectionGWASScenarios.R
SelectionGWASScenariosBinary.R