args <- commandArgs(trailingOnly = TRUE)
#args <- c("../TEMP/PLINKFILES/UKBHapMapSNPsDef22","YearsEducation","../TEMP/YearsEducation.resid.txt","../TEMP/ResultsWeightedYearsEducationChr22.csv","../INPUT/UKBSelectionWeights.csv")

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

#flexiblas_switch(flexiblas_load_backend("/projects/0/Galama/svalten/UKBSelection/CODE/openblas/lib/libopenblas.so"))
flexiblas_set_num_threads(1)
#unpack args
Input<-args[1]
PhenoName<-args[2]
PhenoInput<-args[3]
Output<-args[4]
WeightInput<-args[5]

PhenoResid<-paste0(PhenoName,".resid")
PhenoWeightResid<-paste0(PhenoName,".Wresid")

GenoData <- BEDMatrix(Input, simple_names=TRUE)

SNPlist<-colnames(GenoData)

#Merge with phenotypes, weights here
assert(length(which(duplicated(SNPlist)))==0)
#Conduct GWAS and WGWAS, and placebo GWAS 

Pheno <- fread(PhenoInput,header=TRUE,sep="\ ")

IDCrossWalk <- fread("/projects/0/Galama/andriesm/GWAS_PIPELINE/INPUT/ukb_v3_newbasket.s487395.crosswalk",sep=" ")
#TestAssign <- fread("../TEMP/TestAssign.csv",header=TRUE,sep=",") 
UKBWeights <- fread(WeightInput,header=TRUE,sep=",")

#UKBFamIDs<- fread("../../UKBKin/UKB2_Relatives_withfam.txt", header=TRUE, sep=" ") 
#colnames(UKBFamIDs)<- c("f.eid","FamID")

names(IDCrossWalk)[names(IDCrossWalk) == 'V1'] <- 'IID'
names(IDCrossWalk)[names(IDCrossWalk) == 'V2'] <- 'f.eid'

NROW(Pheno)
Pheno <- join(Pheno,IDCrossWalk, by="f.eid")
NROW(Pheno)
Pheno <- join(Pheno,UKBWeights, by="f.eid")
NROW(Pheno)
#Pheno <- join(Pheno,UKBFamIDs, by="f.eid",match="first")
#NROW(Pheno)

#Pheno$FamID[is.na(Pheno$FamID)]<-0
#Pheno$FamID[Pheno$FamID==0]<-max(Pheno$FamID)+seq(1,NROW(Pheno$FamID[Pheno$FamID==0]),1)

#Better to do this beforehand when cleaning the genetic data.
#Pheno<-Pheno[Pheno$InTraining==0,]
#NROW(Pheno)

Pheno<-Pheno[which(!duplicated(Pheno$f.eid)),]
NROW(Pheno)
Pheno<-Pheno[which(!duplicated(Pheno$IID)),]
NROW(Pheno)
Pheno<-Pheno[which(!is.na(Pheno$IID)),]
Pheno<-Pheno[which(!is.na(Pheno$LassoWeight)),]
NROW(Pheno)

Pheno<-select(Pheno,"IID",all_of(PhenoResid),all_of(PhenoWeightResid),"LassoWeight")
#loop over chunks of 100 SNPs:
Nsnp<-NCOL(GenoData)

print(paste("the number of SNPs test is",Nsnp))

print(paste(detectCores(),"available cores"))

#GWASResults <-data.frame(matrix(NA,nrow=length(SNPlist),ncol=10))
GWASResults <-data.frame(matrix(NA,nrow=0,ncol=6))

#colnames(GWASResults)<-c("SNP","BETA","SE","P","BETA_LassoWT","SE_LassoWT","P_LassoWT","BETA_Placebo","SE_Placebo","P_Placebo")

ptm <- proc.time()
i<-1

for (g in seq(1, Nsnp, 1000)){
 start<-g
 if (g +999 < Nsnp) {end <- g+999} else {end <- Nsnp}
 SNPTempData<-as.data.frame(GenoData[,start:end])
 if (dim(SNPTempData)[2]==1){colnames(SNPTempData)<-colnames(GenoData)[start]}
 
 SNPTempData$IID<-row.names(SNPTempData)
 SNPTempData<-as.data.table(SNPTempData)
 TempData<-merge(SNPTempData,Pheno,by="IID")
 myCluster <- makeCluster(15, # number of cores to use
                          type = "FORK") # type of cluster
 registerDoParallel(myCluster)
 print("cluster made succesfully")
tmp<- foreach(SNP=SNPlist[start:end],.combine=rbind)%dopar%{
   #GWASResults[i,1] <- SNP
   assoc<-lm_robust(formula=formula(paste(PhenoResid,"~",SNP)),se_type="HC0",data=TempData, return_vcov=FALSE)
   #assoc<-lm_robust(formula=formula(paste(PhenoResid,"~",SNP)),se_type="CR0",cluster=FamID,data=TempData, return_vcov=FALSE)
   assoc$coefficients[2]
   assoc$std.error[2]
   assoc$p.value[2]

   assocWT<-lm_robust(formula=formula(paste(PhenoWeightResid,"~",SNP)),se_type="HC0",data=TempData,weights=TempData$LassoWeight, return_vcov=FALSE)
   #assocWT<-lm_robust(formula=formula(paste(PhenoWeightResid,"~",SNP)),se_type="CR0",cluster=FamID,data=TempData,weights=TempData$LassoWeight, return_vcov=FALSE)
   assocWT$coefficients[2]
   assocWT$std.error[2]
   assocWT$p.value[2]
   
  data.frame(assoc$coefficients[2], assoc$std.error[2],assoc$p.value[2],
     assocWT$coefficients[2],assocWT$std.error[2],assocWT$p.value[2])
}
#print(tmp)
GWASResults<-rbind(GWASResults,tmp)
stopCluster(myCluster)
}

proc.time() - ptm


GWASResults$SNP<-rownames(GWASResults)
colnames(GWASResults)<-c("BETA","SE","P","BETA_LassoWT","SE_LassoWT","P_LassoWT","SNP")
GWASResults<-GWASResults[,c("SNP","BETA","SE","P","BETA_LassoWT","SE_LassoWT","P_LassoWT")]

GWASoutput<-paste0("../TEMP/ChromResults/GWAS/",Output)
WGWASoutput<-paste0("../TEMP/ChromResults/WGWAS/",Output)

print(paste("the number of lines found in GWAS output:",NROW(GWASoutput)))
print(paste("the number of lines found in WGWAS output:",NROW(WGWASoutput)))

#Format GWAS summary statistics and export 
BimFile<-read.table(paste0(Input,".bim"), sep ="\t")
colnames(BimFile)<-c("CHR",'SNP',"POS","BP","A1","A2")
BimFile<-BimFile[,c("CHR",'SNP',"BP","A1","A2")]
GWASResults<-join(GWASResults,BimFile,by="SNP")
GWASResults<-GWASResults[,c("SNP","CHR","BP","A1","A2","BETA","SE","P","BETA_LassoWT","SE_LassoWT","P_LassoWT")]

write.table(GWASResults[,c("SNP","CHR","BP","A1","A2","BETA","SE","P")], sep="\t",GWASoutput,row.names=FALSE,quote=FALSE)

WGWASResults<-GWASResults[,c("SNP","CHR","BP","A1","A2","BETA_LassoWT","SE_LassoWT","P_LassoWT")]
names(WGWASResults)<-c("SNP","CHR","BP","A1","A2","BETA","SE","P")

write.table(WGWASResults, WGWASoutput,row.names=FALSE, sep="\t",quote=FALSE)

