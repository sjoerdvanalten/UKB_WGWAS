library(data.table)
args <- commandArgs(trailingOnly = TRUE)

pheno<-as.character(args[1])


data<-fread(paste0("../OUTPUT/Simulations/",pheno,".pheno.glm.linear"),header=TRUE)
data<-as.data.frame(data)
data<-data[c("ID","REF","ALT","OBS_CT","BETA","SE","P")]

data$P<-as.numeric(data$P)

summary(data$P)

NROW(data)
data<-data[!is.na(data$BETA),]
data<-data[!is.na(data$P),]
data$P[data$P==0]<-0.00000000000001
NROW(data)



write.table(data,file=paste0("../TEMP/Simulations/",pheno,"_clean.txt"),sep="\t",col.names=TRUE,row.names=FALSE,quote= FALSE)



 

library(data.table)

  