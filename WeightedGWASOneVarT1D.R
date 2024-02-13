#Import top hits for T1D:

options(bitmapType="cairo")
library(qqman)
library(fastman)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(bigstatsr)
library(estimatr)
library(BEDMatrix)
library(data.table)
library(plyr)

tophits<-read.table("../OUTPUT/TABLES/DifferentSNPsWGWASSigType1Diabetes.csv",sep=",",header=TRUE)
tophitsH<-read.table("../TEMP/TopHitsHausman.csv",sep=",",header=TRUE) #clumped SNPs

tophits<-tophits[tophits$SNP %in% tophitsH$SNP,]

Pheno <- fread("../TEMP/Type1Diabetes.resid.txt",header=TRUE,sep="\ ")
PhenoResid<-"Type1Diabetes.resid"
PhenoWeightResid<-"Type1Diabetes.Wresid"

GWASControls<-readRDS("../TEMP/GWASControls.rds")

IDCrossWalk <- fread("/projects/0/Galama/andriesm/GWAS_PIPELINE/INPUT/ukb_v3_newbasket.s487395.crosswalk",sep=" ")

names(IDCrossWalk)[names(IDCrossWalk) == 'V1'] <- 'IID'
names(IDCrossWalk)[names(IDCrossWalk) == 'V2'] <- 'f.eid'

OrigWeights <- fread("../INPUT/UKBWeightsKFolded.csv",header=TRUE,sep=",")

NROW(Pheno)
Pheno <- join(Pheno,IDCrossWalk, by="f.eid")
Pheno <- join(Pheno,OrigWeights, by="f.eid")
NROW(Pheno)

Pheno<-join(Pheno,GWASControls,by="f.eid")

vars<-c("Education","region","birthyear","sex","carsnoc",
        "HealthSelfReport","Empstat","Tenure","SingleHousehold","Educationregion",
        "HealthSelfReportEducation","EmpstatEducation","HealthSelfReportEmpstat","HealthSelfReportEducationregion","EmpstatEducationregion",
        "EmpstatEducationHealthSelfReport","HealthSelfReportEducationregionEmpstat")
labels<-c("Unweighted","All Variables",
          "Education","Region","Year of Birth","Sex","Number of Cars",
          "Self-reported Health","Household Tenure","Employment Status","Single Household","Education*Region",
          "Self-reported Health * Education","Employment Status * Education","Self-reported Health * Employment Status","Self-reported Health * Education * Region","Empstat * Education * Region",
          "Self-reported Health, Employment Status * Education","Employment Status * Education * Self-reported Health")
color<-c("Unweighted","Weighted",
         "Weighted, selected variables","Weighted, selected variables","Weighted, selected variables","Weighted, selected variables","Weighted, selected variables",
         "Weighted, selected variables","Weighted, selected variables","Weighted, selected variables","Weighted, selected variables","Weighted, selected variables",
         "Weighted, selected variables","Weighted, selected variables","Weighted, selected variables","Weighted, selected variables","Weighted, selected variables",
         "Weighted, selected variables","Weighted, selected variables")

controls<-"Sex+Birthyear+GeneBatch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20"

for (i in 1:NROW(tophits)){
chr<-tophits$CHR[i]
SNP<-tophits$SNP[i]

GenoData <- BEDMatrix(paste0("../TEMP/PLINKFILES/UKBHapMapSNPsDef",chr), simple_names=TRUE)

GenoData<-as.data.frame(GenoData[,SNP])
names(GenoData)<-SNP

GenoData$IID<-row.names(GenoData)
Data<-merge(GenoData,Pheno,by="IID")

robresults <- data.frame(matrix(nrow = 8, ncol = 4)) 


assoc<-lm_robust(formula=formula(paste(PhenoResid,"~",SNP)),se_type="HC0",data=Data, return_vcov=FALSE)
print(summary(assoc))
robresults[1,]<-summary(assoc)$coefficients[SNP,1:4]

assocWT<-lm_robust(formula=formula(paste(PhenoWeightResid,"~",SNP)),se_type="HC0",data=Data,weights=Data$LassoWeight, return_vcov=FALSE)
print(summary(assocWT))
robresults[2,]<-summary(assocWT)$coefficients[SNP,1:4]

assocNoCont<-lm_robust(formula=formula(paste("Type1Diabetes","~",SNP)),se_type="HC0",data=Data, return_vcov=FALSE)
print(summary(assocNoCont))
robresults[3,]<-summary(assocNoCont)$coefficients[SNP,1:4]


assocWTNoCont<-lm_robust(formula=formula(paste("Type1Diabetes","~",SNP)),se_type="HC0",data=Data,weights=Data$LassoWeight, return_vcov=FALSE)
print(summary(assocWTNoCont))
robresults[4,]<-summary(assocWTNoCont)$coefficients[SNP,1:4]

#test robustness against some other specifications:
assoclogit<-glm(formula=formula(paste("Type1Diabetes","~",SNP,"+",controls)),data=Data, family="binomial")
print(summary(assoclogit))
robresults[5,]<-summary(assoclogit)$coefficients[SNP,1:4]


assocWTlogit<-glm(formula=formula(paste("Type1Diabetes","~",SNP,"+",controls)),data=Data,weights=Data$LassoWeight, family="binomial")
print(summary(assocWTlogit))
robresults[6,]<-summary(assocWTlogit)$coefficients[SNP,1:4]

#test robustness against some other specifications:
assoclogitNoCont<-glm(formula=formula(paste("Type1Diabetes","~",SNP)),data=Data, family="binomial")
print(summary(assoclogitNoCont))
robresults[7,]<-summary(assoclogitNoCont)$coefficients[SNP,1:4]

assocWTlogitNoCont<-glm(formula=formula(paste("Type1Diabetes","~",SNP)),data=Data,weights=Data$LassoWeight, family="binomial")
print(summary(assocWTlogitNoCont))
robresults[8,]<-summary(assocWTlogitNoCont)$coefficients[SNP,1:4]

names(robresults)<-c("Estimate","Std. Error","Z-value","P-value")

robresults$type<-c("unweighted LPM","weighted LPM",
                   "unweighted LPM, no controls","weighted LPM, no controls",
                   "unweighted logit", "weighted logit",
                   "unweighted logit, no controls","weighted logit, no controls")

write.table(robresults,paste0("../OUTPUT/TABLES/robresults",SNP,".csv"),col.names=TRUE,row.names=FALSE,quote=FALSE,sep=",")

results<-rbind(summary(assoc)$coefficients[SNP,],summary(assocWT)$coefficients[SNP,])

for (var in vars){
print(var)
UKBWeights <- fread(paste0("../TEMP/OneVarWeights",var,".csv"),header=TRUE,sep=",")
Data <- join(Data,UKBWeights, by="f.eid")
assocNew<-lm_robust(formula=formula(paste(PhenoWeightResid,"~",SNP)),se_type="HC0",data=Data,weights=Weight, return_vcov=FALSE)
print(summary(assocNew))
Data$Weight<-NULL
results<-rbind(results,summary(assocNew)$coefficients[SNP,])
}

results<-cbind(results,labels)
results<-cbind(results,color)

#Next, for weights based on each variable:
results<-as.data.frame(results)
results$labels<-factor(results$labels,levels=labels)
results$color<-factor(results$color)
results$Estimate<-as.numeric(results$Estimate)
results$"CI Upper"<-as.numeric(results$"CI Upper")
results$"CI Lower"<-as.numeric(results$"CI Lower")

ggplot(results, aes(x=labels, y=Estimate,color=color)) + geom_point() + geom_errorbar(aes(ymin=`CI Upper`, ymax=`CI Lower`)) +
  scale_colour_brewer(name="Model",palette="Set1") + theme_bw()  + geom_hline(yintercept=0) + 
  theme(axis.text.x = element_text(angle = 90)) + ylim(-.01,.01) + xlab("Variables used in IP weights")
ggsave(paste0("../OUTPUT/FIGURES/OneVarT1D",SNP,".pdf"))

}


