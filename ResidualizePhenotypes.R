##Date: 03/05/2022
##Author: Sjoerd van Alten
##Goal: A file to residualize a given UKB phenotype for age, sex, recruitment centre, genotype batches, and 10 PCs

args <- commandArgs(trailingOnly = TRUE)

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


form<-paste(varname,"~","Sex+Birthyear+GeneBatch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20")

lm.out <- lm(formula = as.formula(form), data = df)

wlm.out <- lm(formula = as.formula(form), data = df, weights=df$LassoWeight)


PhenoResid<-paste0(varname,".resid")
PhenoWeightResid<-paste0(varname,".Wresid")

df[,PhenoResid]<-df[,varname]-predict(lm.out,newX=df)
df[,PhenoWeightResid]<-df[,varname]-predict(wlm.out,newX=df)

cor(df[,varname],df[,PhenoResid])

cor(df[,varname],df[,PhenoWeightResid])


write.table(df[,c("f.eid", varname, PhenoResid, PhenoWeightResid)], 
            file = paste0("../TEMP/",varname,".resid.txt"),row.names=FALSE,
            append=FALSE, quote=FALSE, sep = "\ ", na = "NA")

