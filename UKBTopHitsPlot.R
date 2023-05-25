args <- commandArgs(trailingOnly = TRUE)

options(bitmapType="cairo")
library(qqman)
library(fastman)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

tophitfile<-as.character(args[1])
#comma seperated list of top hits positions ranges (to use for exclusion in QCtool)

#check: how many independent hits? --> what do we know about these hits? GWAS Catalogue etc? 

#Create localized manhattan plots (+500KB, only weakly correlated (r>=0.1) SNPs)

#merge in GWAS catalog queries as obtained by running the MAF summary stats through a FUMA pipeline:

GWAScat<-read.table("../INPUT/gwas_catalog_v1.0.2-associations_e107_r2022-07-30.tsv",header=TRUE,sep="\t",quote="")

GWAScat<-GWAScat[c("SNP_GENE_IDS","DATE.ADDED.TO.CATALOG","FIRST.AUTHOR","DISEASE.TRAIT","INITIAL.SAMPLE.SIZE","REPLICATION.SAMPLE.SIZE",
                   "SNPS","SNP_ID_CURRENT","P.VALUE","MAPPED_TRAIT")]

names(GWAScat)[names(GWAScat)=="SNPS"]<-"SNP"
names(GWAScat)[names(GWAScat)=="MAPPED_TRAIT"]<-"Trait"
names(GWAScat)[names(GWAScat)=="DATE.ADDED.TO.CATALOG"]<-"Date"
GWAScat$Trait<-gsub("measurement","",GWAScat$Trait)
GWAScat$Trait<-gsub("self reported","",GWAScat$Trait)
GWAScat$Trait<-gsub("body mass index","BMI",GWAScat$Trait)
GWAScat$Trait<-gsub("household income","income",GWAScat$Trait)
GWAScat$Trait<-gsub("sexual intercourse","sex",GWAScat$Trait)
GWAScat$Trait<-gsub("unipolar depression","depression",GWAScat$Trait)
GWAScat$Trait[grep("asthma",GWAScat$Trait)]<-"asthma"
GWAScat$Trait[grep("alcohol",GWAScat$Trait)]<-"alcohol consumption"
GWAScat$Trait[grep("smoking",GWAScat$Trait)]<-"smoking"
GWAScat$Trait[grep("peptic ulcer disease",GWAScat$Trait)]<-"peptic ulcer disease"
GWAScat$Trait<-gsub("(.*),.*","\\1", GWAScat$Trait) #removes elements form string after comma
#Only keep top hits
NROW(GWAScat)
GWAScat[GWAScat$SNP=="rs144942570",]
GWAScat<-GWAScat[as.numeric(GWAScat$"P.VALUE")<=5*10^(-8),]
GWAScat[GWAScat$SNP=="rs144942570",]

NROW(GWAScat)
GWAScat<-GWAScat[rowSums(is.na(GWAScat)) != ncol(GWAScat),]
NROW(GWAScat)


N<-str_extract(GWAScat$"INITIAL.SAMPLE.SIZE","\\(?[0-9,.]+\\)?")
#N<-str_extract_all(GWAScat$"INITIAL.SAMPLE.SIZE","\\(?[0-9,.]+\\)?")

#N<-t(sapply(N,"[", c(3,4,5,18)))
#N<-apply(N,2,gsub,pattern=",",replacement="")
#N<-apply(N,2,as.numeric)
#N<-apply(N,1,sum,na.rm=TRUE)
GWAScat$N<-N
GWAScat$N<-gsub(",","",GWAScat$N)
GWAScat$N<-as.numeric(GWAScat$N)
GWAScat[GWAScat$SNP=="rs144942570",]

tophits<-read.table(tophitfile,sep=",",header=TRUE)

tophits$P<-tophits$P_LassoWT
tophits$BETA<-tophits$BETA_LassoWT
tophits$SE<-tophits$SE_LassoWT
tophits$P_LassoWT<-NULL
tophits$BETA_LassoWT<-NULL
tophits$SE_LassoWT<-NULL
#MAFdata<-read.table(paste0("../TEMP/MAF/WeightedMAFBootstrapTopHits","3",".txt"), sep="\t",header=TRUE)
data<-read.table("../TEMP/ChromResults/WGWAS/HausmanAroundTophitsChr1.csv",header=TRUE)

bimTMP<-read.table("../TEMP/PLINKFILES/UKBAroundTophitsQC1LD0.1.bim")
bimTMP<-bimTMP[c("V1","V2","V4","V5","V6")]
names(bimTMP)<-c("CHR","SNP","BP","A1","A2")
NROW(data)
data<-merge(data,bimTMP)
NROW(data)

LDInfo<-read.table("../TEMP/PLINKFILES/UKBAroundTophitsQC1.ld",header=TRUE)
  
for (chr in c(3,4,5,6,7,11,12,14,15,17,18)){
  #MAFdataTMP<-read.table(paste0("../TEMP/MAF/WeightedMAFBootstrapTopHits",chr,".txt"), sep="\t",header=TRUE)
  dataTMP<-read.table(paste0("../TEMP/ChromResults/WGWAS/HausmanAroundTophitsChr",chr,".csv"),header=TRUE)
  
  bimTMP<-read.table(paste0("../TEMP/PLINKFILES/UKBAroundTophitsQC",chr,"LD0.1.bim"))
  bimTMP<-bimTMP[c("V1","V2","V4","V5","V6")]
  names(bimTMP)<-c("CHR","SNP","BP","A1","A2")
  print(NROW(dataTMP))
  dataTMP<-merge(dataTMP,bimTMP)
  print(NROW(dataTMP))
    
  data<-rbind(data,dataTMP)
    
  LDInfoTMP<-read.table(paste0("../TEMP/PLINKFILES/UKBAroundTophitsQC",chr,".ld"),header=TRUE)
  LDInfo<-rbind(LDInfo,LDInfoTMP)
}
    
hitList<-list()
    
Pmax<-max(-log10(tophits$P))
Pmax<-ceiling(Pmax)
  
for (i in 1:NROW(tophits)){
  hit<-tophits$SNP[i]
  BP<-tophits$BP[i]
  chr<-tophits$CHR[i]
  BPHigh<-BP+500000
  BPLow<-BP-500000
  tmpData<-data[data$BP>=BPLow&data$BP<=BPHigh&data$CHR==chr,]
  LDSNPs<-LDInfo[LDInfo$SNP_A==hit,]
  print(hit)
  print(NROW(tmpData))
  
  names(LDSNPs)[names(LDSNPs)=="SNP_B"]<-"SNP"
  tmpData<-merge(data,LDSNPs[c("SNP","R2")],by="SNP",all.x=FALSE,all.y=FALSE)

  print(NROW(tmpData))
  tmpData$tophit<-tmpData$SNP==hit
  
  NROW(tmpData)
  tmpData<-merge(tmpData,GWAScat,by="SNP",all.x=TRUE,all.y=FALSE)
  
  tmpData$Trait[is.na(tmpData$Trait)]<-""
  #tmpData<-distinct(tmpData,SNP,Trait ,.keep_all=TRUE)
  
  #when there are duplicate triats, keep one in highest LD:
  tmpData<-tmpData[order(tmpData$R2),]
  NROW(tmpData)
  index<-row.names(tmpData)
  print(tmpData)
  tmpData[duplicated(tmpData[,"Trait"], fromLast=T),"Trait"]<-""
  index[which(!index%in%row.names(tmpData))]
  NROW(tmpData)

  print(tmpData[tmpData$Trait!="",])
  
  png(paste0("../OUTPUT/FIGURES/TopHit_LocalManhattan",hit,".png"), width=800, height=800, res=150)
  manhattan(tmpData, chr="CHR", bp="BP", snp="SNP", p="P",highlight=hit, xlim= c(BPLow,BPHigh), suggestiveline=FALSE)
  dev.off()
  
  #build local manhattan plot from scratch:
  #print(tmpData[tmpData$Trait!="",])
  hitList[[i]]<-tmpData
  #png(paste0("../OUTPUT/FIGURES/MAF_WeightedMAF_LocalManhattanAnnotated",hit,".png"), width=1200, height=800, res=150)
  p<-ggplot(tmpData, aes(x=BP, y=-log10(P))) +
    
    # Show all points
    geom_point(aes(color=R2,shape=as.factor(tophit),size=as.factor(tophit)), alpha=1) +
    scale_shape_manual(values=c(16,4),guide="none")+
    scale_size_manual(values=c(1,3),guide="none")+
    {if(NROW(tmpData)>1)scale_color_gradientn(colors=c("green","red","black"),values=c(0,0.5,0.9999999,1),breaks=c(0.11,0.5,0.99,1),labels=c(0.1,0.5,"",1))} +
    
    # Add label using ggrepel to avoid overlapping
    geom_hline(yintercept=-log10(5*10^{-8}), linetype="dashed", color = "blue") + 
    geom_label_repel( data=tmpData, aes(label=Trait), size=3,max.overlaps = Inf, max.time=20, max.iter=100000, force=2, force_pull=0.5) +
    ylim(0,Pmax+1) +
    # Custom the theme:
    theme_bw() + xlab(paste("Chromosome",chr,"position")) + ylab(expression(-log[10](italic(p))))  
    theme( 
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  print(p)
  ggsave(paste0("../OUTPUT/FIGURES/TopHit_LocalManhattanAnnotated",hit,".pdf"),width=1800, height=1200, dpi=288, units="px")
  #dev.off()
}