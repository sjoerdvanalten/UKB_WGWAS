data1<-read.table("../OUTPUT/TABLES/DifferentSNPsWGWASSuggestiveType1Diabetes.csv",header=TRUE,sep=",")
data2<-read.table("../OUTPUT/TABLES/DifferentSNPsWGWASSuggestiveBreastCancer.csv",header=TRUE,sep=",")

ChrListT1D<- unique(data1$CHR)

WGWASClump1<-read.table(paste0("../TEMP/Clumped/","Type1Diabetes","Clumped",ChrListT1D[1],".clumped"),header=TRUE)


for (c in ChrListT1D[2:length(ChrListT1D)]){
  print(c)
  WGWASClumpTemp<-read.table(paste0("../TEMP/Clumped/","Type1Diabetes","Clumped",c,".clumped"),header=TRUE)
  WGWASClump1<-rbind(WGWASClump1,WGWASClumpTemp)
}

ChrListBC<- unique(data2$CHR)
WGWASClump2<-read.table(paste0("../TEMP/Clumped/","BreastCancer","Clumped",ChrListBC[1],".clumped"),header=TRUE)

#independent hits
#data1<-data1[data1$SNP %in% WGWASClump1$SNP,]
#data2<-data2[data2$SNP %in% WGWASClump2$SNP,]

data<-rbind(data1,data2)

data$A1A2<-paste0(data$A1,"/",data$A2)
data$BPRange[data$CHR<10]<-paste0(0,data$CHR[data$CHR<10],":",sapply(data$BP[data$CHR<10]-500000,max,0),"-",data$BP[data$CHR<10]+500000)
data$BPRange[data$CHR>=10]<-paste0(data$CHR[data$CHR>=10],":",sapply(data$BP[data$CHR>=10]-500000,max,0),"-",data$BP[data$CHR>=10]+500000)
data<-data[order(data$CHR),]
#names(tophits)<-c("SNP", "Chromosome", "Location (Bp)", "A1/A2", "MAF (unweighted)", 
#                  "MAF (IP weighted)","P value","Known UKB GWAS hits","Known UKB GWAS hits for loci in high LD")

write.table(as.list(data$BPRange),"../TEMP/TopHitsBPRangeHausman.txt",sep=" ",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(as.list(data$SNP),"../TEMP/TopHitsSNPRangeHausman.txt",sep=" ",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(data,"../TEMP/TopHitsHausman.csv",sep=",",row.names=FALSE,quote=FALSE)

print(unique(data$CHR))
print(sort(unique(data1$CHR)))
print(sort(unique(data2$CHR)))