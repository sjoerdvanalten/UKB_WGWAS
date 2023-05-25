#Goal: Summarize GWAS/WGWAS regressions into a nice table:
library(xtable)

vars<-c("YearsEducation","BMI","SevereObesity","Height","DrinksPerWeek","BreastCancer")
Label<-c("Years of Education","BMI","Severe Obesity","Height","Drinks Per Week","Breast cancer")


table<-c()
for (v in vars){
  tmp<-read.table(paste0("../TEMP/WGWASReg",v,".txt"),header=TRUE,sep="\t")
  table<-rbind(table,tmp)
}

table<-cbind(Label,table)


table$Coefficient<-round(table$Coefficient,digits=3)
table$CILow<-round(table$CILow,digits=3)
table$CIHigh<-round(table$CIHigh,digits=3)

pTemp<-table$P
pTemp[table$P<1e-09]<-paste0("$",substr(as.character(format(table$P[table$P<1e-09]),digits=16),1,4),"\\times 10^{-",substr(as.character(table$P[table$P<1e-09]),nchar(as.character(table$P[table$P<1e-09]))-1
                                                                                                                                       ,nchar(as.character(table$P[table$P<1e-09]))),"}$")
pTemp[table$P>=1e-09&table$P<=1e-03]<-paste0("$",substr(as.character(format(table$P[table$P>=1e-09&table$P<=1e-03]),digits=16),1,4),"\\times 10^{-",substr(as.character(table$P[table$P>=1e-09&table$P<=1e-03]),nchar(as.character(table$P[table$P>=1e-09&table$P<=1e-03]))
                                                                                                                                                                               ,nchar(as.character(table$P[table$P>=1e-09&table$P<=1e-03]))),"}$")


#pTemp[table$P<1e-08]<-"$<10^{-8}$"
#pTemp[table$P>=1e-08&table$P<=1e-05]<-paste0("$",substr(as.character(table$P[table$P>=1e-08&table$P<=1e-04]),1,4),"\\times 10^{-",substr(as.character(table$P_dif[table$P_dif>=1e-08&table$P_dif<=1e-04]),nchar(as.character(table$P_dif[table$P_dif>=1e-08&table$P_dif<=1e-04]))
                                                                                                                                                            # ,nchar(as.character(table$P_dif[table$P_dif>=1e-08&table$P_dif<=1e-04]))),"}$")


pTemp[table$P>1e-03]<-sprintf("%.5f", round(table$P[table$P>1e-03],digits=5))


table$P<-pTemp


table$Coefficient<-paste0("\\textbf{",table$Coefficient,"} [",table$CILow,";",table$CIHigh,"]")

table<-table[c("Label","Coefficient","P","NSNP")]


xTab<-xtable(table,align = "ll|l|l|l")

colnames(xTab) = c("\\textbf{Phenotype}", "\\textbf{Coefficient} [95\\% CI]", "\\textbf{P}","\\textbf{N}")
print(xTab,include.rownames=FALSE,sanitize.colnames.function = paste0,
      sanitize.text.function= paste0,type="latex",file="../OUTPUT/TABLES/WGWASReg.tex",floating=FALSE)
