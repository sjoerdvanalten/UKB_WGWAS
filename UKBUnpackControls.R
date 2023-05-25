#Goal: Unpack GWAS controls: age, sex, recruitment centre, genotype batches, and 20 PCs
#+ 
library(data.table)
df<-fread("/projects/0/Galama/andriesm/GWAS_PIPELINE/OUTPUT/basket2005284.ukb46101_refresh_pheno.csv.gz")
df<-as.data.frame(df)

vars<-c("f.22001.0.0","f.34.0.0","f.22000.0.0",
        "f.22009.0.1","f.22009.0.2","f.22009.0.3","f.22009.0.4", "f.22009.0.5","f.22009.0.6","f.22009.0.7","f.22009.0.8","f.22009.0.9","f.22009.0.10",
            "f.22009.0.11","f.22009.0.12","f.22009.0.13","f.22009.0.14","f.22009.0.15","f.22009.0.16","f.22009.0.17","f.22009.0.18","f.22009.0.19","f.22009.0.20")

varlabels<-c("Sex","Birthyear","GeneBatch",
             "PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10",
             "PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")

df<-df[c("f.eid",vars)]
names(df)<-c("f.eid",varlabels)

df$Sex<-as.factor(df$Sex)
df$Birthyear<-as.factor(df$Birthyear)
df$GeneBatch<-as.factor(df$GeneBatch)
#df$AssessmentCentre<-as.factor(df$AssessmentCentre)

saveRDS(df,"../TEMP/GWASControls.rds")