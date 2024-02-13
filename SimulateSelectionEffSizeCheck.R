c<-10

SNPeffects<-read.table("../TEMP/TestPheno1.effects", header=TRUE)
GWAS<-read.table("../TEMP/Simul_clean.txt", header=TRUE)

summary(abs(SNPeffects$Effect))
CLUMP<-read.table(paste0("../TEMP/Clumped/YearsEducationClumped",10,".clumped"),header=TRUE)
#CLUMP<-read.table("../TEMP/SimulClump.clumped",header=TRUE) 

GWAS<-GWAS[GWAS$ID %in% CLUMP$SNP,]
NROW(GWAS)

TopHits<-SNPeffects[SNPeffects$Var_Explained>=0.0005,]
#TopHits<-SNPeffects
#TopHits<-SNPeffects[order(SNPeffects$Var_Explained),]

names(GWAS)[names(GWAS)=="ID"]<-"Predictor"
GWAS<-GWAS[c("Predictor","BETA")]
TopHits<-merge(TopHits,GWAS)

lm<-lm(Effect~BETA,data=TopHits)
summary(lm)


#TopHits$Effect<-abs(TopHits$Effect)
#TopHits$BETA<-abs(TopHits$BETA)
#lm<-lm(Effect~BETA,data=TopHits)
#summary(lm)