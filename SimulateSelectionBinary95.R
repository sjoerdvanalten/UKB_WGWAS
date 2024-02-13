#Goal: Simulate some data to illustrate the effects of collider bias in the UKB:
library(tidyverse)
library(ggExtra)
library(stargazer)
library(BEDMatrix)
library(data.table)
library(plyr)
library(ggplot2)


options(bitmapType="cairo")

args <- commandArgs(trailingOnly = TRUE)
#args <- c(10,1,2)

chr<-as.numeric(args[1])
caus<-as.character(args[2])
iter<-as.numeric(args[3])

print(chr)
print(caus)
print(iter)

Data<-read.table(paste0("../TEMP/TestPhenoBinary95",caus,".pheno"))
names(Data)<-c("V1","IID","Pheno")

SNPeffects<-read.table(paste0("../TEMP/TestPhenoBinary95",caus,".effects"), header=TRUE)
PGI<-read.table(paste0("../TEMP/SimulPGIBinary95",caus,".profile"),header=TRUE)

Data<-merge(Data,PGI)

Data$SCORE<-(Data$SCORE-mean(Data$SCORE))/sd(Data$SCORE)
summary(lm(formula = Pheno ~ SCORE, data = Data))

SNPeffects<-SNPeffects[order(-SNPeffects$Effect),]
if (as.numeric(caus) <= 10 & as.numeric(caus)>=1){
  SNPinfo<-SNPeffects[1:as.numeric(caus),]
}else{
  SNPinfo<-SNPeffects[1:10,]
}

SNPinfo<-SNPeffects[SNPeffects$Effect==summary(SNPeffects$Effect)["Max."],]

GenoDataFull <- BEDMatrix(paste0("../TEMP/PLINKFILES/UKBHapMapSNPsDef",chr), simple_names=TRUE)

ind<-SNPinfo[1,"Predictor"]

GenoData<-as.data.frame(GenoDataFull[,ind])
names(GenoData)<-ind

GenoData$IID<-row.names(GenoData)

SNPData<-GenoData[,c("IID",ind)]

Data <- join(Data,SNPData,by="IID")

formula<-as.formula(paste0("Pheno~",ind))

#summary(lm(formula,data=Data))

set.seed(3454353)

N<-NROW(Data)

biaslist1<-c()
biaslist2a<-c()
biaslist2b<-c()



betalist<-c(0,0.2,0.4,0.6,0.8,1)
len<-length(betalist)


biasmat1 = matrix(0, nrow = iter, ncol = len)
biasmat2a = matrix(0, nrow = iter, ncol = len)
biasmat2b = matrix(0, nrow = iter, ncol = len)

SEmat1 = matrix(0, nrow = iter, ncol = len)
SEmat2a = matrix(0, nrow = iter, ncol = len)
SEmat2b = matrix(0, nrow = iter, ncol = len)

j<-0
for (beta in betalist){
  j<-j+1
  for (i in 1:iter){
  print(i)
  print(j)
  Data$Error<-rnorm(N, mean = 0, sd = 1)
  Data$S<-beta*Data$Pheno+Data$Error
  Data$S2a<-beta*Data$Pheno+(beta)*Data$SCORE+Data$Error
  Data$S2b<-beta*Data$Pheno-(beta)*Data$SCORE+Data$Error

  #Scenario 1a SNP and y positively related, selection on y positive:
  Data$UKB <- Data$S
  Data$UKB <- as.factor(Data$UKB > quantile(Data$UKB, .95))

  Data$UKB2a <- Data$S2a
  Data$UKB2a <- as.factor(Data$UKB2a > quantile(Data$UKB2a, .95))

  Data$UKB2b <- Data$S2b
  Data$UKB2b <- as.factor(Data$UKB2b > quantile(Data$UKB2b, .95))

  #FillColors<-c("#E69F00", "#56B4E9")

  print(table(Data$Pheno))
  print(table(Data$Pheno[Data$UKB==TRUE]))
  print(table(Data$Pheno[Data$UKB2a==TRUE]))
  print(table(Data$Pheno[Data$UKB2b==TRUE]))
  
  LMTrue<-lm(formula,data=Data)
  LMScen1<-lm(formula,data=Data[Data$UKB==TRUE,])
  LMScen2a<-lm(formula,data=Data[Data$UKB2a==TRUE,])
  LMScen2b<-lm(formula,data=Data[Data$UKB2b==TRUE,])
  
  coef1<-summary(LMScen1)$coefficients[2,1]
  if (NROW(summary(LMScen2a$coefficients)>=2)){
  coef2a<-summary(LMScen2a)$coefficients[2,1]
  SEmat2a[i,j]<-summary(LMScen2a)$coefficients[2,2]
  }else{
    coef2a<-NA
    SEmat2a[i,j]<-NA
  }
  coef2b<-summary(LMScen2b)$coefficients[2,1]
  
  #bias1<-abs(summary(LMScen1)$coefficients[2,1]-truth)/abs(truth)
  #bias2a<-abs(summary(LMScen2a)$coefficients[2,1]-truth)/abs(truth)
  #bias2b<-abs(summary(LMpdwScen2b)$coefficients[2,1]-truth)/abs(truth)

  
print(summary(LMTrue))
print(summary(LMScen1))
print(summary(LMScen2a))
print(summary(LMScen2b))

biasmat1[i,j]<-coef1
biasmat2a[i,j]<-coef2a
biasmat2b[i,j]<-coef2b

SEmat1[i,j]<-summary(LMScen1)$coefficients[2,2]
#
SEmat2b[i,j]<-summary(LMScen2b)$coefficients[2,2]

#Export list with whether individuals are in ``sampled'' UKB

IDList1<-Data[Data$UKB==TRUE,]
Pheno1<-IDList1[,1:3] 
names(Pheno1)<-c("IID","FID","pheno")
Pheno1<-Pheno1[c("FID","IID","pheno")]
Pheno1$pheno<-Pheno1$pheno+2 #add 2 to the phenotype because otherwise plink2 cannot handle the LPM and will switch to logit
write.table(Pheno1,paste0("../TEMP/Simul1Binary95_",i,j,"_",caus,".pheno"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(Pheno1[,1:2],paste0("../TEMP/SimulSelect1Binary95_",i,j,"_",caus,".txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

IDList2a<-Data[Data$UKB2a==TRUE,]
Pheno2a<-IDList2a[,1:3] 
names(Pheno2a)<-c("IID","FID","pheno")
Pheno2a<-Pheno2a[c("FID","IID","pheno")]
Pheno2a$pheno<-Pheno2a$pheno+2 #add 2 to the phenotype because otherwise plink2 cannot handle the LPM and will switch to logit
write.table(Pheno2a,paste0("../TEMP/Simul2aBinary95_",i,j,"_",caus,".pheno"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(Pheno2a[,1:2],paste0("../TEMP/SimulSelect2aBinary95_",i,j,"_",caus,".txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

IDList2b<-Data[Data$UKB2b==TRUE,]
Pheno2b<-IDList2b[,1:3]
names(Pheno2b)<-c("IID","FID","pheno")
Pheno2b<-Pheno2b[c("FID","IID","pheno")]
Pheno2b$pheno<-Pheno2b$pheno+2 #add 2 to the phenotype because otherwise plink2 cannot handle the LPM and will switch to logit


write.table(Pheno2b,paste0("../TEMP/Simul2bBinary95_",i,j,"_",caus,".pheno"),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(Pheno2b[,1:2],paste0("../TEMP/SimulSelect2bBinary95_",i,j,"_",caus,".txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)

}
}


truth<-summary(LMTrue)$coefficients[2,1]
print(truth)

bias1<-abs(colMeans(biasmat1)-truth)/truth
bias2a<-abs(colMeans(biasmat2a)-truth)/truth
bias2b<-abs(colMeans(biasmat2b)-truth)/truth

output<-as.data.frame(cbind(betalist,bias1,bias2a,bias2b))

biasmat1<-as.data.frame(biasmat1)
biasmat2a<-as.data.frame(biasmat2a)
biasmat2b<-as.data.frame(biasmat2b)

biasmat1$truth<-truth


write.table(biasmat1,paste0("../TEMP/Simul1ResultsBinary95",caus,".txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(biasmat2a,paste0("../TEMP/Simul1Results2aBinary95",caus,".txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(biasmat2b,paste0("../TEMP/Simul1Results2bBinary95",caus,".txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)


write.table(SEmat1,paste0("../TEMP/Simul1SEsBinary95",caus,".txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(SEmat2a,paste0("../TEMP/Simul1SEs2aBinary95",caus,".txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
write.table(SEmat2b,paste0("../TEMP/Simul1SEs2bBinary95",caus,".txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)


#png(paste0("../OUTPUT/FIGURES/Simulbias.png"), width=2880, height=1800, res=288)
#ggplot(data=output,aes(x=betalist,y=bias1)) +
#  geom_line(linetype = "dashed",color='1')+
#  geom_point() + 
#  geom_line(aes(x=betalist,y=bias2a),linetype = "solid",color='2a')+ 
#  geom_point(aes(x=betalist,y=bias2a)) + 
#  geom_line(aes(x=betalist,y=bias2b),linetype = "dotted",color='2b')+ 
#  geom_point(aes(x=betalist,y=bias2b)) + scale_linetype_manual(name='Scenario',
#                                            breaks=c('1', '2a', '2b'),
#                                            values=c('1'='dashed', '2a'='solid', '2b'='dotted'))
#dev.off()
#+theme_bw()

bias1 
bias2a
bias2b
#True std. errors vs. average empirical std. error

print("scen1")
apply(biasmat1,2,sd)
colMeans(SEmat1)

print("scen2a")
apply(biasmat2a,2,sd)
colMeans(SEmat2a)

print("scen2b")
apply(biasmat2b,2,sd)
colMeans(SEmat2b)




