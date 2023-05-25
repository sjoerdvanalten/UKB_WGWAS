rm(list=ls())

library(plyr)
library(dplyr)
library(stargazer)
library(Hmisc)
library(reporttools)
library(pROC)
library(lmtest)
library(sandwich)
library(stats)
library(purrr)
library(testit)
library(data.table)
library(xtable)
library(estimatr)

options(bitmapType="cairo")

WeightedMean<-function(x, w, na.rm=TRUE,CI=FALSE){
  if (length(x)!=length(w)){
    stop("x and weights are not of same length")
  }
  if (na.rm==TRUE){ 
    w <- w[!is.na(x)]
    x <- x[!is.na(x)]
  }
  m<-sum(x*w)/sum(w)
  if (CI==TRUE){
  se=sqrt(sum(w*(x-m)^2)/sum(w))
  CILow<-m-(se/sqrt(NROW(x)))*1.96
  CIHigh<-m+(se/sqrt(NROW(x)))*1.96
  return(cbind(m,CILow,CIHigh))
  }else{return(m)}
  }


WeightedSE<-function(x, w, na.rm=TRUE){
  if (length(x)!=length(w)){
    stop("x and weights are not of same length")
  }
  if (na.rm==TRUE){ 
    w <- w[!is.na(x)]
    x <- x[!is.na(x)]
  }
  m<-sum(x*w)/sum(w)
  se=sqrt(sum(w*(x-m)^2)/sum(w))
  return(se)
}


setwd("/projects/0/Galama/UKB_WGWAS/CODE")

##here

Pheno<- c("actModVig","bmi","childrenEverMothered","height",
    "t1d","ageFirstBirth","cad","highBloodPressure",
    "neuroticismScore","t2d","ageParents90th","cancerBreast","depress",
    "householdIncome","nonCancerIllness","totChol","ageSmoke","cancer",
    "depressScore","insomniaFrequent","obesitySevere","alzheimer","cancerProstate",        
    "dpw","lifeSatisfaction","positiveAffect","wellBeingSpectrum","anxiety",        
    "cataract","educYears","loneliness","risk","worryFeeling",
    "arthritis","cesSmoke","healthRating","maxcpd","smokeInit",
    "asthma","childrenEverFathered","hearingDifficulty","medsTaken","stroke")

PhenoLabels<-c("Physical activity","BMI","Children mothered","Height",
               "Diabetes - Type 1","Age at first birth (female)","Coronary heart disease","High blood pressure",
               "Neuroticism","Diabetes - Type 2","Age parents ($>$90th pct.)","Breast cancer","Depression",
               "Household income","No. of illnesses (not cancer)","Cholesterol","Age started smoking (former smokers)","Cancer",
               "Depression Score","Insomnia","Severe obesity","Alzheimer","Prostate cancer",
               "Drinks per week","Life satisfaction","General happiness","Well-being","Anxiety",
               "Cateract","Years of education","Loneliness","Risk","Worry",
               "Arthritis","Smoking cessation","Health rating","Max. cigarettes per day","Smoking initiation",
               "Asthma","Children ever fathered","Hearing difficulty","Meds taken","Stroke")


categories<-c("Health behavior","Anthropometric","Health behavior","Anthropometric",
              "Health","Health behavior","Health","Health",
              "Personal trait","Health","Personal trait","Health","Mental health",
              "Personal trait","Health","Health","Health behavior","Health",
              "Mental health","Health","Anthropometric","Mental health","Health",
              "Health behavior","Mental health","Mental health","Mental health","Mental health",
              "Health","Personal trait","Mental health","Personal traits","Mental health",
              "Health","Health behavior","Health","Health behavior","Health behavior",
              "Health","Health behavior","Health","Health","Health")


#Load all phenotypes
Phenofolder<-"../INPUT/PHENOTYPES/"
UKB<-read.table(paste0(Phenofolder,Pheno[1],"_pheno.PREPARED.txt"),header=TRUE)
temp<-read.table(paste0(Phenofolder,"actModVig","_pheno.PREPARED.txt"),header=TRUE)
for (var in Pheno[2:length(Pheno)]){
  temp<-read.table(paste0(Phenofolder,var,"_pheno.PREPARED.txt"),header=TRUE)
  UKB<-merge(UKB,temp,all.x=TRUE,all.y=TRUE)
}


temp<-read.table(paste0("../TEMP/BreastCancer.Female.resid.txt"),header=TRUE)

UKB$BreastCancer<-NULL

crosswalk<-read.table("/projects/0/Galama/andriesm/GWAS_PIPELINE/INPUT/ukb_v3_newbasket.s487395.crosswalk")
names(crosswalk)<-c("IID","f.eid")


temp<-merge(temp,crosswalk)


temp<-temp[c("BreastCancer","IID")]

str(temp)
str(UKB)


UKB<-merge(UKB,temp,all.x=TRUE,all.y=TRUE)

NROW(UKB)

head(UKB)

WinsorLassoWeights<-read.table("../INPUT/UKBWeightsKFolded.csv",header=TRUE,sep=",")
IDCrossWalk <- fread("/projects/0/Galama/andriesm/GWAS_PIPELINE/INPUT/ukb_v3_newbasket.s487395.crosswalk",sep=" ")

names(IDCrossWalk)<-c('IID','f.eid')
UKB <- join(UKB,IDCrossWalk, by="IID",match="first")
UKB <- join(UKB,WinsorLassoWeights, by="f.eid",match="first")
N1<-NROW(UKB)
UKB<-UKB[!is.na(UKB$LassoWeight),]
N2<-NROW(UKB)

paste("dropped",N1-N2,"because they had no estimated weights")


N2/N1
g<-c()
for (v in Pheno){
  print(v)
  x<-t.test(UKB[,v])$conf.int[1:2]
  y<-t.test(UKB[,v])$estimate[1]
  m<-c(x,y)
  l<-WeightedMean(UKB[,v],UKB$LassoWeight,na.rm=TRUE,CI=TRUE)
  
  SD<-sd(UKB[,v],na.rm=TRUE)
  SD_W<-WeightedSE(UKB[,v],UKB$LassoWeight)
  N<-NROW(UKB[,v][!is.na(UKB[,v])])
  t<-c(m,SD,l,SD_W,N)
  g<-rbind(g,t)
  }

g<-cbind(PhenoLabels,g)
g<-cbind(Pheno,g)
g<-cbind(categories,g)
g<-as.data.frame(g)
colnames(g)<-c("Category","Var","Varlabels","CILow","CIHigh","Mean","SD","WMean","WCILow","WCIHigh","WSD","N")
g<-g[order(g$Category),]
g$Change<-abs((as.numeric(g$WMean)-as.numeric(g$Mean))/as.numeric(g$Mean))


g<-g[g$Varlabels %in% c("BMI","Height","Severe obesity","Diabetes - Type 1","Breast cancer","Health rating",
                        "Physical activity","Age at first birth (female)","Drinks per week","Years of education"), ] #comment out to get all phenotype. This is just the selection for the preprint.

cat<-unique(g$Category)
catindex<-which(!duplicated(g$Category))-1

g$Mean<-as.numeric(g$Mean)
g$CILow<-as.numeric(g$CILow)
g$CIHigh<-as.numeric(g$CIHigh)
g$WMean<-as.numeric(g$WMean)
g$WCILow<-as.numeric(g$WCILow)
g$WCIHigh<-as.numeric(g$WCIHigh)
g$SD<-as.numeric(g$SD)
g$WSD<-as.numeric(g$WSD)

g$Mean<-round(g$Mean,digits=3)
g$WMean<-round(g$WMean,digits=3)
g$CILow<-round(g$CILow,digits=3)
g$CIHigh<-round(g$CIHigh,digits=3)
g$WCILow<-round(g$WCILow,digits=3)
g$WCIHigh<-round(g$WCIHigh,digits=3)
g$Change<-round(g$Change,digits=3)
g$SD<-round(g$SD,digits=3)
g$WSD<-round(g$WSD,digits=3)

g$Mean<-paste0("\\textbf{",g$Mean,"} [",g$CILow,";",g$CIHigh,"]")
g$WeightedMean<-paste0("\\textbf{",g$WMean,"} [",g$WCILow,"; ",g$WCIHigh,"]")
g$Change<-paste0(g$Change*100,"\\%")

g<-g[c("Varlabels","Mean","SD","WeightedMean","WSD","N","Change")]

# \\% change after weighting each mean calculated as the absolute value of the weighted mean, minus the unweighted mean, divided by the unweighted mean.
SumTab<-xtable(g,align = "ll|l|l|l|l|l|l", caption="\\textbf{Mean of various phenotypes in the UKB before and after weighting using IPWs:} 95\\% confidence intervals around each mean included.",label="tab:pheno")
addtorow<-list()
addtorow$pos<-as.list(catindex)
cat<-paste0("\\hline \\textit{",cat,"} \\\\ \\hline ")
addtorow$command<-as.vector(cat)

colnames(SumTab) = c("\\textbf{Variable}", "\\textbf{Mean [95\\% CI]}", "\\textbf{SD}", "\\textbf{Weighted Mean [95\\% CI]}", "\\textbf{Weighted SD}","\\textbf{N}","\\textbf{\\% Mean change}")
print(SumTab,include.rownames=FALSE,add.to.row=addtorow,sanitize.colnames.function = paste0,
      sanitize.text.function= paste0,type="latex",file="../OUTPUT/TABLES/WeightedPhenotypesDescriptives.tex",floating=FALSE,size="scriptsize",tabular.environment="longtable")

print(mean(as.numeric(g$N)))

write.table(g,file="../OUTPUT/TABLES/WeightedPhenotypesDescriptives.csv",sep=",",row.names=FALSE,quote=FALSE)


#Again, for breast cancer and sex stratified

SexData<-read.table("../TEMP/sex.txt",header=TRUE)[c("f.eid","Sex")]

BUKB<-UKB[c("f.eid","cancerBreast","LassoWeight")]
BUKB<-merge(BUKB,SexData)

  g<-c()
  
  x<-t.test(BUKB[,"cancerBreast"])$conf.int[1:2]
  y<-t.test(BUKB[,"cancerBreast"])$estimate[1]
  m<-c(x,y)
  
  l<-WeightedMean(BUKB[,"cancerBreast"],BUKB$LassoWeight,na.rm=TRUE,CI=TRUE)
  SD<-sd(BUKB[,"cancerBreast"],na.rm=TRUE)
  SD_W<-WeightedSE(BUKB[,"cancerBreast"],UKB$LassoWeight)
  N<-NROW(BUKB[,"cancerBreast"][!is.na(BUKB[,"cancerBreast"])])
  t<-c(m,SD,l,SD_W,N)
  g<-rbind(g,t)

  x<-t.test(BUKB[BUKB$Sex==1,"cancerBreast"])$conf.int[1:2]
  y<-t.test(BUKB[BUKB$Sex==1,"cancerBreast"])$estimate[1]
  m<-c(x,y)
  
  l<-WeightedMean(BUKB[BUKB$Sex==1,"cancerBreast"],BUKB$LassoWeight[BUKB$Sex==1],na.rm=TRUE,CI=TRUE)
  SD<-sd(BUKB[BUKB$Sex==1,"cancerBreast"],na.rm=TRUE)
  SD_W<-WeightedSE(BUKB[BUKB$Sex==1,"cancerBreast"],UKB$LassoWeight[BUKB$Sex==1])
  N<-NROW(BUKB[BUKB$Sex==1,"cancerBreast"][!is.na(BUKB[BUKB$Sex==1,"cancerBreast"])])
  t<-c(m,SD,l,SD_W,N)
  g<-rbind(g,t)
  
  x<-t.test(BUKB[BUKB$Sex==0,"cancerBreast"])$conf.int[1:2]
  y<-t.test(BUKB[BUKB$Sex==0,"cancerBreast"])$estimate[1]
  m<-c(x,y)
  
  l<-WeightedMean(BUKB[BUKB$Sex==0,"cancerBreast"],BUKB$LassoWeight[BUKB$Sex==0],na.rm=TRUE,CI=TRUE)
  SD<-sd(BUKB[BUKB$Sex==0,"cancerBreast"],na.rm=TRUE)
  SD_W<-WeightedSE(BUKB[BUKB$Sex==0,"cancerBreast"],UKB$LassoWeight[BUKB$Sex==0])
  N<-NROW(BUKB[BUKB$Sex==0,"cancerBreast"][!is.na(BUKB[BUKB$Sex==0,"cancerBreast"])])
  t<-c(m,SD,l,SD_W,N)
  g<-rbind(g,t)
  
  
g<-cbind(c("Total","Male","Female"),g)
g<-as.data.frame(g)
colnames(g)<-c("label","CILow","CIHigh","Mean","SD","WMean","WCILow","WCIHigh","WSD","N")

g$Mean<-as.numeric(g$Mean)
g$CILow<-as.numeric(g$CILow)
g$CIHigh<-as.numeric(g$CIHigh)
g$WMean<-as.numeric(g$WMean)
g$WCILow<-as.numeric(g$WCILow)
g$WCIHigh<-as.numeric(g$WCIHigh)
g$SD<-as.numeric(g$SD)
g$WSD<-as.numeric(g$WSD)

g$Mean<-format(round(g$Mean,digits=4),scientific=FALSE)
g$WMean<-format(round(g$WMean,digits=4),scientific=FALSE)
g$CILow<-format(round(g$CILow,digits=4),scientific=FALSE)
g$CIHigh<-format(round(g$CIHigh,digits=4),scientific=FALSE)
g$WCILow<-format(round(g$WCILow,digits=4),scientific=FALSE)
g$WCIHigh<-format(round(g$WCIHigh,digits=4),scientific=FALSE)
g$SD<-format(round(g$SD,digits=4),scientific=FALSE)
g$WSD<-format(round(g$WSD,digits=4),scientific=FALSE)

g$Mean<-paste0("\\textbf{",g$Mean,"} [",g$CILow,";",g$CIHigh,"]")
g$WeightedMean<-paste0("\\textbf{",g$WMean,"} [",g$WCILow,"; ",g$WCIHigh,"]")

g<-g[c("label","Mean","SD","WeightedMean","WSD","N")]
g<-xtable(g)

print.xtable(g,file="../OUTPUT/TABLES/BreastCancerStratifiedSummary.tex",floating=FALSE,include.rownames=FALSE,
             sanitize.text.function=function(x){x})

