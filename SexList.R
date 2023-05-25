GWASControls<-readRDS("../TEMP/GWASControls.rds")

#Make column names Sex.Resid and Sex.WResid, simply because the WeightedGWAS script can't handle it otherwise (these are not residualized)
varname<-"Sex"


crosswalk<-read.table("/projects/0/Galama/andriesm/GWAS_PIPELINE/INPUT/ukb_v3_newbasket.s487395.crosswalk")
WeightData<-read.table("../INPUT/UKBWeightsKFolded.csv",header=TRUE,sep=",")

names(crosswalk)<-c("FID","f.eid")
GWASControls<-merge(GWASControls,crosswalk)
df<-GWASControls
df<-merge(df,WeightData)

df$Sex<-as.numeric(df$Sex)-1
summary(df$Sex)

form<-paste(varname,"~Birthyear+GeneBatch+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20")

lm.out <- lm(formula = as.formula(form), data = df)

wlm.out <- lm(formula = as.formula(form), data = df, weights=df$LassoWeight)


PhenoResid<-paste0(varname,".resid")
PhenoWeightResid<-paste0(varname,".Wresid")


df[,PhenoResid]<-df[,varname]-predict(lm.out,newX=df)
df[,PhenoWeightResid]<-df[,varname]-predict(wlm.out,newX=df)

cor(df[,varname],df[,PhenoResid])

cor(df[,varname],df[,PhenoWeightResid])

#female=1
write.table(df[,c("f.eid", varname, PhenoResid, PhenoWeightResid)],
            file = paste0("../TEMP/sex.txt"),row.names=FALSE,
            append=FALSE, quote=FALSE, sep = "\ ", na = "NA")

