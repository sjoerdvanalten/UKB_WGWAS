options(bitmapType="cairo")
args <- commandArgs(trailingOnly = TRUE)

library(qqman)
library(fastman)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stats)
library(readxl)
library(data.table)

args<-c("../TEMP/GWASonTheWeights","../TEMP/GWASOnTheWeightsClumped","GWASOnTheWeights")
GWASPre<-as.character(args[1]) #../TEMP/GWASonTheWeights
ClumpPre<-as.character(args[2]) #../TEMP/GWASOnTheWeightsClumped
outfix<-as.character(args[3]) #GWASOnTheWeights

results<-read.table(paste0(GWASPre,"1.P1.assoc.linear"),header=TRUE)
bim<-read.table("../TEMP/PLINKFILES/UKBHapMapSNPsDef1.bim")
names(bim)<-c("CHR","SNP","POS","BP","A1","A2")
clump<-read.table(paste0(ClumpPre,"1.txt.clumped"),header=TRUE)


for (chr in c(rep(2:22))){
  print(chr)
  resultsTMP<-read.table(paste0(GWASPre,chr,".P1.assoc.linear"), header=TRUE)
  bimTMP<-read.table(paste0("../TEMP/PLINKFILES/UKBHapMapSNPsDef",chr,".bim"))
  clumpTMP<-read.table(paste0(ClumpPre,chr,".txt.clumped"),header=TRUE)
  names(bimTMP)<-c("CHR","SNP","POS","BP","A1","A2")
  results<-rbind(results,resultsTMP)
  bim<-rbind(bim,bimTMP)
  clump<-rbind(clump,clumpTMP)
}

results<-results[results$TEST!="INTERCEPT",]
#Attach correct A2 to results file
bim$A2_GWAS<-NA
bim$A2_GWAS[bim$A1==results$A1]<-bim$A2[bim$A1==results$A1]
bim$A2_GWAS[bim$A1!=results$A1]<-bim$A1[bim$A1!=results$A1]
bim$A2<-bim$A2_GWAS

results<-merge(results,bim[c("SNP","A2")])

results<-results[c("SNP","CHR","BP","A1","A2","TEST","NMISS","BETA","STAT","P")]

write.table(results,file=paste0("../TEMP/",outfix,"Full.tab"),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

#p<-ggplot(tmpData, aes(x=BP, y=-log10(PBootstrap))) +
#  
#  # Show all points
#  geom_point(aes(color=R2,shape=as.factor(tophit),size=as.factor(tophit)), alpha=1) +
#  scale_shape_manual(values=c(16,4),guide="none")+
#  scale_size_manual(values=c(1,3),guide="none")+
#  scale_color_gradientn(colors=c("green","red","black"),values=c(0,0.5,0.9999999,1),breaks=c(0.11,0.5,0.99,1),labels=c(0.1,0.5,"",1)) +
  
  # Add label using ggrepel to avoid overlapping
#  geom_hline(yintercept=-log10(5*10^{-8}), linetype="dashed", color = "blue") + 
#  geom_label_repel( data=tmpData, aes(label=Trait), size=3,max.overlaps = Inf, max.time=20, max.iter=100000, force=2, force_pull=0.5) +
#  ylim(0,Pmax+1) +
  # Custom the theme:
#  theme_bw() + xlab(paste("Chromosome",chr,"position")) + ylab(expression(-log[10](italic(p))))  
#theme( 
#  panel.border = element_blank(),
#  panel.grid.major.x = element_blank(),
#  panel.grid.minor.x = element_blank()
#)
#print(p)
#ggsave(paste0("../OUTPUT/FIGURES/MAF_WeightedMAF_LocalManhattanAnnotated",hit,".pdf"),width=2400, height=1600, dpi=288, units="px")
#dev.off()
#}

pdf(paste0("../OUTPUT/FIGURES/",outfix,"_QQ.pdf"),width=7,height=7)
fastqq(results$P)
dev.off()

nrow(results)
resultsIndep <- merge(results,clump[c("SNP")],by="SNP")
nrow(resultsIndep)

resultsIndep[resultsIndep$P<5*10^(-8),]

TopListMerge<-resultsIndep[resultsIndep$P<5*10^(-8),]
TopListMerge$TopHit<-"yes"
TopListMerge<-TopListMerge[c("SNP","TopHit")]

NROW(results)
results<-merge(results,TopListMerge,all=TRUE)
NROW(results)

results$TopHit[is.na(results$TopHit)]<-"no"
#Make manhattan plot
#pdf(paste0("../OUTPUT/FIGURES/",outfix,"_Manhattan.pdf"))
#manhattan(results, chr="CHR", bp="BP", snp="SNP", p="P",annotateTop=TRUE,suggestiveline=FALSE)
#dev.off() 
don <- results %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(results, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)


axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

don<-don[don$P<0.1,]
Pmax<-max(-log10(don$P))

ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("black", "grey"), 22 )) +
  geom_hline(yintercept=-log10(5*10^(-8)), linetype="dashed", 
             color = "red") + 
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  
  #annotate top hits
  geom_label_repel( data=subset(don, TopHit=="yes"), aes(label=SNP), size=3,max.overlaps = Inf, max.time=20, max.iter=100000, force=2, force_pull=0.5) +
  #scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  ylim(0,Pmax+1) +
  # Custom the theme:
  theme_bw() +  xlab(paste("Base pair position")) + ylab(expression(-log[10](italic(p))))  +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
ggsave(paste0("../OUTPUT/FIGURES/",outfix,"_Manhattan.pdf"),width=18,height=12,units="cm")


#pleiotropy analysis: 
pleio<-read_excel("../INPUT/Wattanabe2019Supplementary.xlsx", sheet = "ST 4",skip=1)
#pleio$`#domains`[pleio$`#domains`>5]<-5
names(pleio)[names(pleio)=="#domains"]<-"pleiotraits"
pleio<-pleio[c("chr","start","end","pleiotraits","#loci")]
#merge data with nearest starting value in pleio data set. Next, check if BP is within the range of the pleio locus. If not, set `pleiotraits`==0
DT<-resultsIndep
DT$startMerge<-DT$BP
#DT<-DT[,c("CHR","start")]
DT <- data.table(DT, key = c("CHR","startMerge"))

tm<-pleio[c("chr","start","end","pleiotraits")]
tm$startMerge<-tm$start
names(tm)[names(tm)=="chr"]<-"CHR"
tm <- data.table(tm, key = key(DT))
NROW(DT)
DTest<-tm[DT, roll='nearest']#merge on chromosome and start, but based on nearest for start
NROW(DTest)
DTest$startMerge<-NULL

DTest$pleiotraits[DTest$BP>=DTest$start&DTest$BP<=DTest$end]<-NA
table(DTest$pleiotraits)

mean(DTest$pleiotraits, na.rm=TRUE)
mean(DTest$pleiotraits[DTest$P<=5*10^(-8)], na.rm=TRUE)

DTest$Top<-as.numeric(DTest$P<=5*10^(-8))
NROW(DTest)
DTest<-DTest[DTest$P>=5*10^(-5)|DTest$P<5*10^(-8),]
NROW(DTest)

table(DTest$Top)

t.test(pleiotraits ~ Top, data = DTest)
 

DTest$pleiotraits<-as.factor(DTest$pleiotraits)
obs <- table(DTest$pleiotraits[DTest$P<=5*10^(-8)])/NROW(DTest[DTest$P<=5*10^(-8)*!is.na(DTest$pleiotraits),])
exp <- table(DTest$pleiotraits)/length(DTest$pleiotraits[!is.na(DTest$pleiotraits)])
# Main test
chisq.test(obs, p=exp,correct=TRUE)$p.value

#merge in base-pair position from .bim file 
summary(resultsIndep$P)
paste("the number of independent, genome-wide significant hits is", NROW(resultsIndep$P[resultsIndep$P<5*10^(-8)]))
tophits<-resultsIndep[resultsIndep$P<5*10^(-8),]

suggestivehits<-resultsIndep[resultsIndep$P<5*10^(-5),]
paste("the number of independent, suggestive hits is", NROW(suggestivehits))


tophits$A1A2<-paste0(tophits$A1,"/",tophits$A2)
tophits$BPRange[tophits$CHR<10]<-paste0(0,tophits$CHR[tophits$CHR<10],":",sapply(tophits$BP[tophits$CHR<10]-500000,max,0),"-",tophits$BP[tophits$CHR<10]+500000)
tophits$BPRange[tophits$CHR>=10]<-paste0(tophits$CHR[tophits$CHR>=10],":",sapply(tophits$BP[tophits$CHR>=10]-500000,max,0),"-",tophits$BP[tophits$CHR>=10]+500000)
tophits<-tophits[order(tophits$CHR),]
#names(tophits)<-c("SNP", "Chromosome", "Location (Bp)", "A1/A2", "MAF (unweighted)", 
#                  "MAF (IP weighted)","P value","Known UKB GWAS hits","Known UKB GWAS hits for loci in high LD")
write.table(tophits,paste0("../OUTPUT/TABLES/TopHits",outfix,".csv"),sep=",",row.names=FALSE,quote=FALSE)

suggestivehits$BPRange[suggestivehits$CHR<10]<-paste0(0,suggestivehits$CHR[suggestivehits$CHR<10],":",sapply(suggestivehits$BP[suggestivehits$CHR<10]-500000,max,0),"-",suggestivehits$BP[suggestivehits$CHR<10]+500000)
suggestivehits$BPRange[suggestivehits$CHR>=10]<-paste0(suggestivehits$CHR[suggestivehits$CHR>=10],":",sapply(suggestivehits$BP[suggestivehits$CHR>=10]-500000,max,0),"-",suggestivehits$BP[suggestivehits$CHR>=10]+500000)
suggestivehits<-suggestivehits[order(suggestivehits$CHR),]
suggestivehits$A1A2<-paste0(suggestivehits$A1,"/",suggestivehits$A2)

write.table(suggestivehits,paste0("../OUTPUT/TABLES/SuggestiveHits",outfix,".csv"),sep=",",row.names=FALSE,quote=FALSE)

write.table(as.list(tophits$BPRange),"../TEMP/TopHitsBPRangeGWASOnTheWeights.txt",sep=" ",row.names=FALSE,quote=FALSE,col.names=FALSE)
write.table(as.list(tophits$SNP),"../TEMP/TopHitsSNPRangeGWASOnTheWeights.txt",sep=" ",row.names=FALSE,quote=FALSE,col.names=FALSE)

table(tophits$CHR)

#Median effect size here, possibly correcting for winner's curse. 
summary(resultsIndep$BETA)

png("../OUTPUT/FIGURES/GWASOnTheWeightsEffectSizes.png", units="px", width=800, height=800, res=150)
hist(resultsIndep$BETA)
dev.off()
