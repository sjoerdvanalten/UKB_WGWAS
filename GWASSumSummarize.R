#Goal: Summarize GWAS/WGWAS results into a nice table:
library(xtable)

vars<-c("YearsEducation","BMI","Height","DrinksPerWeek","SevereObesity","AgeFirstBirth","BreastCancer","PhysicalActivity","Type1Diabetes","HealthRating")
Label<-c("Years of Education","BMI","Height","Drinks per Week","Severe Obesity","Age at First Birth","Breast cancer",
          "Physical activity","Type 1 Diabetes","Self-rated health")
table<-c()
for (v in vars){
  tmp<-read.table(paste0("../TEMP/WGWASSum",v,".txt"),header=TRUE,sep="\t")
  table<-rbind(table,tmp)
}

table<-cbind(Label,table)

LD<-read.table("../TEMP/LDScoreResults.txt",header=TRUE,sep="\t")

LD$Phenotype[LD$Phenotype=="BreastCancerFemale"]<-"BreastCancer"

LD$P<-(1-pnorm(abs((LD$rG-1)/LD$rGSE)))/2

LD$rG<-round(LD$rG,3)

LD$rG[LD$P<0.005]<-paste0("$",LD$rG[LD$P<0.005],"^{*}")
LD$rG[LD$P<0.005]<-paste0(LD$rG[LD$P<0.005]," (",LD$rGSE[LD$P<0.005],")$")
LD$rG[LD$P>=0.005]<-paste0(LD$rG[LD$P>=0.005]," (",LD$rGSE[LD$P>=0.005],")")
LD<-LD[c("Phenotype","rG")]

table<-merge(table,LD)

table$SEIncrease<-paste0(round(table$SEIncrease*100,1),"\\%")

names(table)[names(table)=="rG"]<-"rG_GWAS_WGWAS"
table<-table[c("Label","rG_GWAS_WGWAS","NeffGWAS","NeffIPWGWAS","SEIncrease","yPGWAS5e.08","yPWGWAS5e.08","nUnique","nUniqueSig")]

tableX<-rbind(c(rep(" ",NCOL(table))),table)

xTab<-xtable(tableX,align = "ll|l|l|l|l|l|l|l|l", caption="\\textbf{Caption here:} describe p-values")

colnames(xTab) = c("\\multirow{2}{*}{\\textbf{Phenotype}}", "\\multirow{2}{*}{\\textbf{$r(\\hat{\\beta}_{GWAS},\\hat{\\beta}_{WGWAS})$}}", "\\multirow{2}{*}{$N^{GWAS}_{eff}$}", "\\multirow{2}{*}{$N^{WGWAS}_{eff}$}",
                   "\\multirow{2}{*}{Increase \\\\ S.E.s}","\\multirow{2}{*}{\\shortstack[l]{sig. hits \\\\ GWAS}}","\\multirow{2}{*}{\\shortstack[l]{sig. hits \\\\ WGWAS}}","\\multirow{2}{*}{\\shortstack[l]{$sig. hits \\\\ WGWAS, unique$}}",
                   "\\multirow{2}{*}{\\shortstack[l]{New loci}}")
print(xTab,include.rownames=FALSE,sanitize.colnames.function = paste0,
      sanitize.text.function= paste0,type="latex",file="../OUTPUT/TABLES/GWASSumSummarize.tex",hline.after=c(-1,1,nrow(xTab)),floating=FALSE)

summary(table$NeffGWAS)
summary(table$NeffIPWGWAS)

1-mean(table$NeffIPWGWAS)/mean(table$NeffGWAS)