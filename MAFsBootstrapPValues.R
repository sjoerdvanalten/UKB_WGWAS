#Goal: Calculate weighted MAFs, including Bootstrap for the variance of the difference between weighted and unweighted MAFs.
library(flexiblas)

args <- commandArgs(trailingOnly = TRUE)

library(flexiblas)
library(data.table)
library(plyr)
library(testit)
library(svMisc)
library(estimatr)
library(tidyverse)
library(doParallel)
library(foreach)
library(bigstatsr)
library(BEDMatrix)
library(modi)
library(dplyr)

flexiblas_set_num_threads(1)

Input<-as.character(args[1])
Chr<-as.character(args[2])
cores<-args[3]

set.seed(9832)
WeightInput<-"../INPUT/UKBWeightsKFolded.csv"
print(Input)
GenoData <- BEDMatrix(Input, simple_names=TRUE)

#Uncomment to sample 100 SNPs only (for debugging purposes):
#selection<-sample(colnames(GenoData),100)
#selection<-colnames(GenoData)
#GenoData<-GenoData[,selection]

SNPlist<-colnames(GenoData)
Nsnp<-NCOL(GenoData)

UKBWeights <- fread(WeightInput,header=TRUE,sep=",")
IDCrossWalk <- fread("/projects/0/Galama/andriesm/GWAS_PIPELINE/INPUT/ukb_v3_newbasket.s487395.crosswalk",sep=" ")

names(IDCrossWalk)[names(IDCrossWalk) == 'V1'] <- 'IID'
names(IDCrossWalk)[names(IDCrossWalk) == 'V2'] <- 'f.eid'

#loop over chunks of 1000 SNPs:
Nsnp<-NCOL(GenoData)

print(Nsnp)

MAF <- matrix(nrow=0,ncol=3)

rm(GenoData) #need to load GenoData within each core
#bootstrap the difference of MAF-MAF_w
MAF<-bootstrapMAF<-function(SNPName,data,iterations=100){
  set.seed(3650)
  BootstrapDif<-c(1:iterations)
  for(j in 1:iterations){
    TempData<-slice_sample(data,n=NROW(data),replace=TRUE)
    TempMAF<-mean(TempData[,SNPName],na.rm=TRUE)/2 
    TempMAFW<-weighted.mean(TempData[,SNPName],w=TempData$LassoWeight,na.rm=TRUE)/2
    BootstrapDif[j]<-TempMAF-TempMAFW 
  }
  return(var(BootstrapDif))
}

print(paste(detectCores(),"available cores"))

ptm <- proc.time()

i<-1
myCluster <- makeCluster(cores, # number of cores to use
                         type = "FORK") # type of cluster
registerDoParallel(myCluster)
print("cluster made succesfully")

#Read the genetic data in in chunks of 100 SNPs, calculated, maf, weighted maf, and the variance for their
#+difference, for each SNP:
MAF<-foreach(g=seq(1, Nsnp, 100),.combine=rbind)%dopar%{
  GenoData <- BEDMatrix(Input, simple_names=TRUE)
  start<-g
  if (g +99 < Nsnp) {end <- g+99} else {end <- Nsnp}
  SNPTempData<-as.data.frame(GenoData[,start:end])
  if (dim(SNPTempData)[2]==1){colnames(SNPTempData)<-colnames(GenoData)[start]}
  tempN<-end-start+1
  SNPTempData$IID<-row.names(SNPTempData)
  SNPTempData <- join(SNPTempData,IDCrossWalk, by="IID")
  SNPTempData <- join(SNPTempData,UKBWeights, by="f.eid")
  SNPTempData <- SNPTempData[!is.na(SNPTempData$LassoWeight),]
  maf<-apply(SNPTempData[,1:tempN],2,mean,na.rm=TRUE)
  maf<-round(maf/2,digits=4)  
  maf_weighted<-apply(SNPTempData[,1:tempN],2,weighted.mean,w=SNPTempData$LassoWeight,na.rm=TRUE)
  maf_weighted<-round(maf_weighted/2,digits=4)  
  varBootstrap<-sapply(colnames(SNPTempData[,1:tempN]),bootstrapMAF,data=SNPTempData)
  cbind(maf,maf_weighted,varBootstrap)   
  }


proc.time() - ptm

NROW(MAF)

MAF<-as.data.frame(MAF)
MAF$ZBootstrap<-(MAF$maf-MAF$maf_weighted)/(sqrt(MAF$varBootstrap)) #calculate Z-stat
MAF$PBootstrap<-2*pnorm(-abs(MAF$ZBootstrap)) #Calculate P-value
MAF$SNP<-row.names(MAF)
MAF<-MAF[c("SNP","maf","maf_weighted","varBootstrap","ZBootstrap","PBootstrap")]

write.table(MAF,file=paste0("../TEMP/MAF/WeightedMAFBootstrap",Chr,".txt"), sep="\t", row.names=FALSE,quote=FALSE)

