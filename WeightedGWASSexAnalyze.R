args <- commandArgs(trailingOnly = TRUE)

#args<-c("YearsEducationChr","YearsEducation","../../svalten/GWAS_SUMMARY/EduYears_Main.txt.hapmap","MarkerName","A1","A2","BETA","SE","../TEMP/GWASSum/EA2Chr","../TEMP/YearsEducation.resid.txt","YearsEducation")
#args<-c("BMIChr","BMI","../../svalten/GWAS_SUMMARY/SNP_gwas_mc_merge_nogc.tbl.uniq.hapmap","SNP","A1","A2","BETA","SE","../TEMP/GWASSum/BMI1Chr","../TEMP/BMI.resid.txt","BMI")
#args<-c("HeightChr","Height","../../svalten/GWAS_SUMMARY/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.hapmap","MarkerName","Allele1","Allele2","BETA","SE","../TEMP/GWASSum/Height1Chr")
#args<-c("AgeFirstBirthChr","AgeFirstBirth","../../svalten/GWAS_SUMMARY/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.hapmap","MarkerName","Allele1","Allele2","BETA","SE","../TEMP/GWASSum/Height1Chr")
#args<-c("DepressionChr","Depression","NA","NA","NA","NA","NA","NA","NA","../TEMP/Depression.resid.txt")
#args<-c("Sex","Sex","../../svalten/GWAS_SUMMARY/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.hapmap","MarkerName","Allele1","Allele2","BETA","SE","../TEMP/GWASSum/Height1Chr","../TEMP/sex.txt")
#args<-c("BreastCancerChr","BreastCancer","../TEMP/oncoarray_bcac_public_release_oct17_clean.txt.hapmap","SNP","a0","a1","BETA","SE","../TEMP/GWASSum/BreastCancerChr","../TEMP/BreastCancer.resid.txt","BreastCancer")
#args<-c("SexChr","Sex",NA,NA,NA,NA,NA,NA,NA,"../TEMP/sex.txt","Sex")

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

library(qqman)
library(fastman)
library(testit)
library(plyr)
library(lmtest)
library(ggplot2)
library(AER)
library(data.table)

options(bitmapType="cairo")

#unpack args
Input<-args[1]
PhenoName<-args[2]
OriginalGWAS<-args[3]
MarkerName<-args[4]
A1Name<-args[5]
A2Name<-args[6]
BetaName<-args[7]
SEName<-args[8]
Clumpfile<-args[9]
Phenofile<-as.character(args[10])
ClumpPrefix<-as.character(args[11])

#Results<-read.table(paste0("../TEMP/ChromResults/GWAS/",Input,1,".tab"),sep="\t",header=TRUE)
Results<-read.table(paste0("../TEMP/ChromResults/GWAS/",Input,1,".tab"),sep="\t",header=TRUE)

WResults<-read.table(paste0("../TEMP/ChromResults/WGWAS/",Input,1,".tab"),sep="\t",header=TRUE)

MAF<-read.table("../TEMP/MAF/WeightedMAF1.txt", header=TRUE)
Pheno<-read.table(Phenofile,header=TRUE)
WeightInput<-"../INPUT/UKBWeightsKFolded.csv"

IDCrossWalk <- fread("/projects/0/Galama/andriesm/GWAS_PIPELINE/INPUT/ukb_v3_newbasket.s487395.crosswalk",sep=" ")
UKBWeights <- fread(WeightInput,header=TRUE,sep=",")
Pheno<-join(Pheno,IDCrossWalk)
Pheno<-join(Pheno,UKBWeights)

s_sq<-var(Pheno[,paste0(PhenoName,".resid")],na.rm=TRUE)
s_wsq<-WeightedSE(Pheno[,paste0(PhenoName,".resid")],Pheno$LassoWeight,na.rm=TRUE)^2

GWASClump<-read.table(paste0("../TEMP/Clumped/",ClumpPrefix,"Clumped",1,".clumped"),header=TRUE)
WGWASClump<-read.table(paste0("../TEMP/Clumped/",ClumpPrefix,"Clumped",1,".clumped"),header=TRUE)

for (c in 2:22){
  print(c)
  temp<-read.table(paste0("../TEMP/ChromResults/GWAS/",Input,c,".tab"),sep="\t",header=TRUE)
  Results<-rbind(Results,temp)
  Wtemp<-read.table(paste0("../TEMP/ChromResults/WGWAS/",Input,c,".tab"),sep="\t",header=TRUE)
  WResults<-rbind(WResults,Wtemp)
  Mtemp<-read.table(paste0("../TEMP/MAF/WeightedMAF",c,".txt"),header=TRUE)
  MAF<-rbind(MAF,Mtemp)
  GWASClumpTemp<-read.table(paste0("../TEMP/Clumped/",ClumpPrefix,"Clumped",c,".clumped"),header=TRUE)
  WGWASClumpTemp<-read.table(paste0("../TEMP/Clumped/",ClumpPrefix,"Clumped",c,".clumped"),header=TRUE)
  GWASClump<-rbind(GWASClump,GWASClumpTemp)
  WGWASClump<-rbind(WGWASClump,WGWASClumpTemp)
  }

NROW(Results)

names(WResults)<-c("SNP","CHR","BP","A1","A2","BETA_LassoWT","SE_LassoWT","P_LassoWT")
Results<-merge(Results,WResults)
head(Results)

NROW(Results)

Results<-merge(Results,MAF)
WResults<-merge(WResults,MAF)

#effective sample size, estimate using formula in within-sibship GWAS (Howe et al.,):
Results$N_eff<-s_sq/(Results$SE^2*(2*Results$MAF*(1-Results$MAF)))
WResults$N_eff<-s_wsq/(Results$SE_LassoWT^2*(2*Results$MAF_weighted*(1-Results$MAF_weighted)))

names(WResults)[names(WResults)=="BETA_LassoWT"]<-"BETA"
names(WResults)[names(WResults)=="SE_LassoWT"]<-"SE"
names(WResults)[names(WResults)=="P_LassoWT"]<-"P"

N_effW<-mean(WResults$N_eff)
N_effGWAS<-mean(Results$N_eff)

Results$Z<-Results$BETA/Results$SE
WResults$Z<-WResults$BETA/WResults$SE

write.table(Results[c("SNP","CHR","BP","A1","A2","BETA","SE","P","MAF","N_eff","Z")],paste0("../OUTPUT/GWAS/",ClumpPrefix,".txt"),quote=FALSE,row.names=FALSE)
write.table(WResults[c("SNP","CHR","BP","A1","A2","BETA","SE","P","MAF","N_eff","Z")],paste0("../OUTPUT/WGWAS/",ClumpPrefix,".txt"),quote=FALSE,row.names=FALSE)

WResults<-NULL

PValuesReport<-function(Pvec){
  result<-list(NA,NA)
  names(result)<-c("5e-08","5e-06")
  count<-1
  for (threshold in c(5*10^(-8),5*10^(-6))){
    SigNo<-length(Pvec[Pvec<threshold])
    result[count]<-SigNo
    print(paste("Number of genomewide significant p-values at",threshold,"equals",SigNo))
    count<-count+1
  }
  return(result)
  #return(data.frame(Data$Predictor[Data$Wald_P<5*10^(-8)])) 
}

#Keep independent hits only.
GWASPTestVec<-Results$P[Results$SNP %in% GWASClump$SNP]
NROW(GWASPTestVec)
WGWASPTestVec<-Results$P_LassoWT[Results$SNP %in% WGWASClump$SNP]
NROW(WGWASPTestVec)


yPGWAS<-PValuesReport(GWASPTestVec)
yPWGWAS<-PValuesReport(WGWASPTestVec)

#Check by how much standard errors increase as a result of our weighting procedure:
Results$SE_WTincrease<-(Results$SE_LassoWT-Results$SE)/(Results$SE)

SEIncrease<-mean(Results$SE_WTincrease)

xCorr<-cor.test(Results$BETA,Results$BETA_LassoWT)$estimate
xCorrCI<-cor.test(Results$BETA,Results$BETA_LassoWT)$conf.int[1:2]
print(xCorr)

SumResultsOut<-unlist(c(PhenoName,xCorr,round(N_effGWAS),round(N_effW),SEIncrease,yPGWAS,yPWGWAS))
SumResultsOut<-as.list(SumResultsOut)
#Make p-values for the difference in the SNP:

SumResultsOut<-as.data.frame(SumResultsOut)

names(SumResultsOut)<-c("Phenotype","Correlation","NeffGWAS","NeffIPWGWAS","SEIncrease","yPGWAS5e-08","yPGWAS5e-06","yPWGWAS5e-08","yPWGWAS5e-06")
write.table(SumResultsOut,paste0("../TEMP/WGWASSum",PhenoName,".txt"),row.names=FALSE,quote=FALSE,sep="\t")

print("Hausman P-values")

Results$VarDifHausman<-Results$SE_LassoWT^2-Results$SE^2
Results$PHausman<-2*pnorm(-abs((Results$BETA-Results$BETA_LassoWT)/sqrt(Results$VarDifHausman)))

summary(Results$PHausman)

summary(Results$BETA)
summary(Results$BETA_LassoWT)

print(head(Results[Results$PHausman<5*10^(-8),]))

pdf(paste0("../OUTPUT/FIGURES/QQDif",PhenoName,".pdf"))
fastqq(Results$PHausman)
dev.off()

if (OriginalGWAS!="NA"){

#Compare to the original GWAS:
OrigResults <- read.table(OriginalGWAS,header=TRUE,sep="\t")
head(OrigResults)
OrigResults <- OrigResults[c(MarkerName,BetaName,SEName,A1Name,A2Name)]
names(OrigResults)<-c('SNP','BETA_orig','SE_orig','A1_orig','A2_orig')
head(OrigResults)
#Results<-Results[Results$BETA>-0.25,]

nrow(OrigResults)
OrigResults<-OrigResults[OrigResults$A1_orig=="A"|OrigResults$A1_orig=="G"|OrigResults$A1_orig=="C"|OrigResults$A1_orig=="T",]
nrow(OrigResults)
OrigResults<-OrigResults[OrigResults$A2_orig=="A"|OrigResults$A2_orig=="G"|OrigResults$A2_orig=="C"|OrigResults$A2_orig=="T",]
nrow(OrigResults)

Results <- merge(Results, OrigResults, by="SNP")

#only clumped SNPs

ClumpResults <- read.table(paste0(Clumpfile,1,".clumped"),header=TRUE)
for (c in 2:22){
  if (file.exists(paste0(Clumpfile,c,".clumped"))){
  Ctemp<-read.table(paste0(Clumpfile,c,".clumped"),header=TRUE)
  ClumpResults<-rbind(ClumpResults,Ctemp)
  }
}

Results <- merge(Results,ClumpResults[c("SNP")],by="SNP")

#NullResults<-read.table(paste0("../OUTPUT/ResultsLasso",PhenoName,"FullAllClumped.txt"),sep="\t",header=TRUE)
#NROW(NullResults)
#NullResults<-NullResults[NullResults$P>0.5,]
#NROW(NullResults)
#NullResults<-NullResults[NullResults$P_LassoWT>0.5,]
#NROW(NullResults)  
#Re<-cor(NullResults$BETA,NullResults$BETA_LassoWT)

#PValuesReport(Results$P_Placebo)

xSE<-mean(Results$SE_WTincrease,na.rm=TRUE)

Results$BETA_orig[Results$A1==Results$A2_orig&Results$A1!=Results$A1_orig]<--Results$BETA_orig[Results$A1==Results$A2_orig&Results$A1!=Results$A1_orig]
Results$A1_orig[Results$A1==Results$A2_orig&Results$A1!=Results$A1_orig]<-Results$A1[Results$A1==Results$A2_orig&Results$A1!=Results$A1_orig]
Results$A2_orig[Results$A1_orig==Results$A2_orig]<-Results$A2[Results$A1_orig==Results$A2_orig]

NROW(Results)
Results<-Results[Results$A1==Results$A1_orig,]
NROW(Results)
assert(Results$A1==Results$A1_orig)
assert(Results$A2==Results$A2_orig)

cor(Results$BETA,Results$BETA_orig)
cor(Results$BETA_LassoWT,Results$BETA_orig)

Results$Z_orig<-Results$BETA_orig/Results$SE_orig
summary(Results$Z_orig)
N10e05<-NROW(Results$Z_orig[abs(Results$Z_orig)>4.417])

LassoLM<-lm(BETA_LassoWT~BETA-1,data=Results)
LassoLM_196<-lm(BETA_LassoWT~BETA-1,data=Results[abs(Results$Z_orig)>1.96,])
LassoLM_1e05<-lm(BETA_LassoWT~BETA-1,data=Results[abs(Results$Z_orig)>4.417,])
LassoLM_1e07<-lm(BETA_LassoWT~BETA-1,data=Results[abs(Results$Z_orig)>5.326724,])


summary(LassoLM)
summary(LassoLM_196)
summary(LassoLM_1e05)
summary(LassoLM_1e07)

NModel<-NROW(Results[abs(Results$Z_orig)>4.417,])

CILassoLM<-coefci(LassoLM_1e05)

xP<-linearHypothesis(LassoLM_1e05, c("BETA=1"))$"Pr(>F)"[2]

xCoef<-LassoLM_1e05$coef
xCI<-CILassoLM

RegResults<-as.list(c(PhenoName,xCoef,xCI,xP,NModel))
print(RegResults)
RegResults<-as.data.frame(RegResults)
print(RegResults)
names(RegResults)<-c("Phenotype","Coefficient","CILow","CIHigh","P","NSNP")
write.table(RegResults,paste0("../TEMP/WGWASReg",PhenoName,".txt"),row.names=FALSE,quote=FALSE,sep="\t")

p<-ggplot(Results[abs(Results$Z_orig)>4.417,], aes(x=BETA, y=BETA_LassoWT)) + 
  geom_errorbar(aes(ymin = BETA_LassoWT-1.96*SE_LassoWT, ymax = BETA_LassoWT+1.96*SE_LassoWT), colour="#827d7d", alpha=0.5)+ 
  geom_errorbar(aes(xmin = BETA-1.96*SE, xmax = BETA+1.96*SE), colour="#827d7d", alpha=0.5)  +
  geom_point()  + 
  geom_smooth(method=lm, se=TRUE, formula=y~x-1) + geom_abline(intercept=0,slope=1,linetype="dashed") +
  annotate(x=mean(Results$BETA), y=max(Results[abs(Results$Z_orig)>4.417,]$BETA_LassoWT+1.6*Results[abs(Results$Z_orig)>4.417,]$SE_LassoWT), 
           label=paste("slope ==", round(LassoLM_1e05$coefficients[1],3)," (CI ==",round(CILassoLM[1,1], 3),",",round(CILassoLM[1,2],3),")"), 
           geom="text", parse=TRUE,size=6) + ylab(expression(paste(hat(beta)[IPWGWAS], ' (corrected for volunteering)'))) +
  annotate(x=mean(Results$BETA), y=0.8*max(Results[abs(Results$Z_orig)>4.417,]$BETA_LassoWT+1.6*Results[abs(Results$Z_orig)>4.417,]$SE_LassoWT), 
           label=paste(NModel,"SNPs"), geom="text",size=6) + 
  ylab(expression(paste(hat(beta)[IPWGWAS], ' (corrected for volunteering)')))+ xlab(expression(hat(beta)[OLS])) + theme_bw() +
  ylim(min(Results[abs(Results$Z_orig)>4.417,]$BETA_LassoWT-2*Results[abs(Results$Z_orig)>4.417,]$SE_LassoWT),max(Results[abs(Results$Z_orig)>4.417,]$BETA_LassoWT+2*Results[abs(Results$Z_orig)>4.417,]$SE_LassoWT)) +
  xlim(min(Results[abs(Results$Z_orig)>4.417,]$BETA_LassoWT-2*Results[abs(Results$Z_orig)>4.417,]$SE_LassoWT),max(Results[abs(Results$Z_orig)>4.417,]$BETA_LassoWT+2*Results[abs(Results$Z_orig)>4.417,]$SE_LassoWT))+ 
  theme(axis.text = element_text(size=20), axis.title=element_text(size=20)) 
ggsave(plot=p,height=6,width=6,dpi=200, filename=paste0("../OUTPUT/FIGURES/Scatter",PhenoName,"TopHits.pdf"))
#dev.off()

}


#xOut<-c(PhenoName,xCorr,xCorrCI,xSE,xCoef,xCI,xP,xNSNPs,xNConcord)
#names(xOut)<-c("PhenoName","Corr","CorrCILow","CorrCIHigh","SEIncrease","RegCoef","RegCILow","RegCIHigh","RegP","NSNP","Concordance")

#saveRDS(xOut,paste0("../TEMP/",PhenoName,"WGWASSum.rds"))

#yOut<-cbind(yPGWAS,yPWGWAS)

#saveRDS(yOut,paste0("../TEMP/",PhenoName,"WGWASPValues.rds"))
