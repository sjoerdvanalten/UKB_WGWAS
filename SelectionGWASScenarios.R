#Goal: Simulate some data to illustrate the effects of collider bias in the UKB:
library(tidyverse)
library(ggExtra)
library(stargazer)

setwd("C:/Users/sjoerd/Dropbox/UKB_WGWAS/CODE")

set.seed(3454353)

N<-1000000
MAF<-0.4
SNP<-rbinom(N,2,MAF) #simulate SNP data with MAF of 0.4
table(SNP)
Error<-rnorm(N, mean = 0, sd = 1)

Y<-0.04*SNP+Error

#Scenario 1a SNP and y positively related, selection on y positive:
Data<-data.frame(SNP,Y,Error)
Data$UKB <- Data$Y
Data$UKB <- as.factor(Data$UKB > quantile(Data$UKB, .95))
Data$UKB_under <- as.factor(Data$Y < quantile(Data$Y, .05))

FillColors<-c("#E69F00", "#56B4E9")

LMTrue<-lm(Y~SNP,data=Data)
LMScen1<-lm(Y~SNP,data=Data[Data$UKB==TRUE,])
LMScen1_under<-lm(Y~SNP,data=Data[Data$UKB_under==TRUE,])

MAFTrue<-t.test(Data$SNP/2)

MAFScen1<-t.test(Data$SNP[Data$UKB==TRUE]/2)

#Scenario 2a, SNP and y positivelY related, selection on SNP and y positive:
Data<-data.frame(SNP,Y,Error)
Data$UKBt <- 0.5*Data$Y+0.04*Data$SNP
Data$UKB <- as.factor(Data$UKBt > quantile(Data$UKBt, .95))
Data$UKB_under <- as.factor(Data$UKBt < quantile(Data$UKBt, .05))


Data$UKBt<-NULL

MAFScen2a<-t.test(Data$SNP[Data$UKB==TRUE]/2)

LMScen2a<-lm(Y~SNP,data=Data[Data$UKB==TRUE,])
LMScen2a_under<-lm(Y~SNP,data=Data[Data$UKB_under==TRUE,])

t.test(Data$SNP[Data$UKB==TRUE]/2)

#Scenario 2b, y positively, SNP negatively related, selection on SNP and y positive:
Data<-data.frame(SNP,Y,Error)
Data$UKB <- 0.5*Data$Y-0.04*Data$SNP
Data$UKB <- as.factor(Data$UKB > quantile(Data$UKB, .95))

mean(Data$SNP)/2

mean(Data$SNP[Data$UKB==TRUE])/2

LMScen2b<-lm(Y~SNP,data=Data[Data$UKB==TRUE,])


MAFScen2b<-t.test(Data$SNP[Data$UKB==TRUE]/2)



addline<-list(c("MAF (95\\% CI)", 
                paste0(round(MAFTrue$estimate, digits = 3),
                      "[",
                      round(MAFTrue$conf.int[1], digits = 4),":",
                      round(MAFTrue$conf.int[2], digits = 4), 
                      "]"),
                paste0(round(MAFScen1$estimate, digits = 3), 
                      "[",
                      round(MAFScen1$conf.int[1], digits = 4),":",
                      round(MAFScen1$conf.int[2], digits = 4), 
                      "]"),
                paste0(round(MAFScen2a$estimate, digits = 3),
                      "[",
                      round(MAFScen2a$conf.int[1], digits = 4),":",
                      round(MAFScen2a$conf.int[2], digits = 4), 
                      "]"),
                paste0(round(MAFScen2b$estimate, digits = 3),
                      "[",
                      round(MAFScen2b$conf.int[1], digits = 4),":",
                      round(MAFScen2b$conf.int[2], digits = 4), 
                      "]")))
names(addline)<-NULL

stargazer(LMTrue,LMScen1,LMScen2a,LMScen2b,type="html",out="../OUTPUT/TABLES/SimulationResults.doc",
          column.labels=c("Population","Scenario 1","Scenario 2a","Scenario 2b"),df=FALSE,keep.stat=c("n","rsq"),
          omit.table.layout = "n",star.cutoffs = NA)


stargazer(LMTrue,LMScen1,LMScen2a,LMScen2b,out="../OUTPUT/TABLES/SimulationResults.tex",
          column.labels=c("Population","Scenario 1","Scenario 2a","Scenario 2b"),df=FALSE,
          dep.var.caption="Phenotype",dep.var.labels.include = FALSE,keep.stat=c("n","rsq"),
          add.lines = addline,
          omit.table.layout = "n",star.cutoffs = NA, float=FALSE)


