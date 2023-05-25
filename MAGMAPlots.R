library(ggplot2)
library(stringr)

setwd("C:/Users/sjoerd/Dropbox/UKB_WGWAS")

for (pheno in c("AFB","BMI","BreastCancer_female","DrinksPerWeek","HealthRating",
                "Height","PhysicalActivity","SevereObesity","Type1Diabetes","YearsEducation")){
#unweighted results

tab1<-read.table(paste0("OUTPUT/MAGMA/",pheno,"_unweighted/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out"),header=TRUE)

tab1$FULL_NAME<-str_replace_all(tab1$FULL_NAME,"_"," ")
tab1<-tab1[order(tab1$P),]
tab1$P<--log10(tab1$P)
tab1$FULL_NAME<-factor(tab1$FULL_NAME)
tab1<-tab1[1:24,]
Plim1<-max(tab1$P)+1

tab<-read.table(paste0("OUTPUT/MAGMA/",pheno,"_weighted/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out"),header=TRUE)

tab$FULL_NAME<-str_replace_all(tab$FULL_NAME,"_"," ")
tab<-tab[order(tab$P),]
tab$P<--log10(tab$P)
tab$FULL_NAME<-factor(tab$FULL_NAME)
tab<-tab[1:24,]
Plim2<-max(tab$P)+1

Plim<-max(Plim1,Plim2)


tab1$sig<-tab1$P>-log10(0.05/54)
p<-ggplot(data=tab1, aes(x=reorder(FULL_NAME, -P),y=P,fill=sig))+geom_bar(stat="identity",show.legend=FALSE)+
  geom_hline(yintercept=-log10(0.05/54),linetype="dashed") + scale_fill_manual(values=c("#2227b6","#990000")) +
  ylab("-log 10 P-value") + xlab("")+theme_bw() + scale_y_continuous(limits = c(0, Plim))+
  theme(axis.text.x = element_text(angle = 70, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename=paste0("OUTPUT/FIGURES/MAGMA/",pheno,".png"),width=6, height=4)

tab$sig<-tab$P>-log10(0.05/54)
p<-ggplot(data=tab, aes(x=reorder(FULL_NAME, -P),y=P,fill=sig))+geom_bar(stat="identity",show.legend=FALSE)+
  geom_hline(yintercept=-log10(0.05/54),linetype="dashed") + scale_fill_manual(values=c("#2227b6","#990000")) +
  ylab("-log 10 P-value") + xlab("")+theme_bw() + scale_y_continuous(limits = c(0, Plim))+
  theme(axis.text.x = element_text(angle = 70, hjust=1),panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename=paste0("OUTPUT/FIGURES/MAGMA/",pheno,"_weighted.png"),width=6, height=4)
}