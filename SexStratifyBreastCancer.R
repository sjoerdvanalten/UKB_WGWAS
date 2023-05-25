
args <- commandArgs(trailingOnly = TRUE)
args<-c("../INPUT/PHENOTYPES/cancerBreast_pheno.PREPARED.txt","BreastCancer")
library(dplyr)
library(plyr)
library(data.table)

input<-as.character(args[1])
varname<-as.character(args[2])

print(varname)

GWASControls<-readRDS("../TEMP/GWASControls.rds")
crosswalk<-read.table("/projects/0/Galama/andriesm/GWAS_PIPELINE/INPUT/ukb_v3_newbasket.s487395.crosswalk")
pheno<-read.table(input,header=TRUE)
WeightData<-read.table("../INPUT/UKBWeightsKFolded.csv",header=TRUE,sep=",")

names(crosswalk)<-c("FID","f.eid")
GWASControls<-merge(GWASControls,crosswalk)

names(pheno)[3]<-varname
df<-merge(pheno,GWASControls)
df<-merge(df,WeightData)

dfMale<-df[df$Sex==1,]
dfFemale<-df[df$Sex==0,]
rm(df)
form<-paste(varname,"~","Birthyear+GeneBatch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20")

lmMale.out <- lm(formula = as.formula(form), data = dfMale)

wlmMale.out <- lm(formula = as.formula(form), data = dfMale, weights=dfMale$LassoWeight)

PhenoResid<-paste0(varname,".resid")
PhenoWeightResid<-paste0(varname,".Wresid")

dfMale[,PhenoResid]<-dfMale[,varname]-predict(lmMale.out,newX=dfMale)
dfMale[,PhenoWeightResid]<-dfMale[,varname]-predict(wlmMale.out,newX=dfMale)

cor(dfMale[,varname],dfMale[,PhenoResid])

cor(dfMale[,varname],dfMale[,PhenoWeightResid])

write.table(dfMale[,c("f.eid", varname, PhenoResid, PhenoWeightResid)], 
            file = paste0("../TEMP/",varname,".Male.resid.txt"),row.names=FALSE,
            append=FALSE, quote=FALSE, sep = "\ ", na = "NA")


lmFemale.out <- lm(formula = as.formula(form), data = dfFemale)

wlmFemale.out <- lm(formula = as.formula(form), data = dfFemale, weights=dfFemale$LassoWeight)

dfFemale[,PhenoResid]<-dfFemale[,varname]-predict(lmFemale.out,newX=dfFemale)
dfFemale[,PhenoWeightResid]<-dfFemale[,varname]-predict(wlmFemale.out,newX=dfFemale)

cor(dfFemale[,varname],dfFemale[,PhenoResid])

cor(dfFemale[,varname],dfFemale[,PhenoWeightResid])

write.table(dfFemale[,c("f.eid", varname, PhenoResid, PhenoWeightResid)], 
            file = paste0("../TEMP/",varname,".Female.resid.txt"),row.names=FALSE,
            append=FALSE, quote=FALSE, sep = "\ ", na = "NA")

