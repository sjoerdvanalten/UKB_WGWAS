Weights<-read.table("../INPUT/UKBWeightsKFolded.csv",header=TRUE,sep=",")

Crosswalk<-read.table("/projects/0/Galama/andriesm/GWAS_PIPELINE/INPUT/ukb_v3_newbasket.s487395.crosswalk")
names(Crosswalk)<-c("IID","f.eid")

Weights<-merge(Weights,Crosswalk)

#order according to .fam file:
fam<-read.table("../TEMP/PLINKFILES/UKBHapMapSNPsDef1.fam")
fam<-fam[c("V1","V2")]
names(fam)<-c("FID","IID")

NROW(Weights)
Weights<-merge(Weights,fam)
NROW(Weights)

pheno<-Weights[c("FID","IID","LassoWeight")]
write.table(pheno,file="../TEMP/UKBWeightsPheno.pheno",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

sample<-read.table("../TEMP/samplesbolt.sample")
sample<-sample[3:NROW(sample),]

#Residualize from controls:
GWASControls<-readRDS("../TEMP/GWASControls.rds")

Weights<-merge(Weights,GWASControls)

names(sample)<-c("FID","IID")

Weights<-merge(Weights,pheno)


form<-"LassoWeight~Sex+Birthyear+GeneBatch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20"

lm.out <- lm(formula = as.formula(form), data = Weights)
Weights[,"LassoWeightResid"]<-Weights[,"LassoWeight"]-predict(lm.out,newX=Weights)

cor(Weights[,"LassoWeight"],Weights[,"LassoWeightResid"])

Weights$FID<-Weights$IID

pheno<-Weights[c("FID","IID","LassoWeightResid")]


head(pheno)

write.table(pheno,file="../TEMP/UKBWeightsPhenoResid.pheno",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
write.table(pheno,file="../TEMP/UKBWeightsPhenoResid.pheno_bolt",col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")