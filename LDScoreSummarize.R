#Goal: Summarize LDScore Regresion results into a nice table:
library(xtable)

PDifference<-function(h2,Wh2,h2se,Wh2se,Intercept){
  Z<-(h2-Wh2)/sqrt(h2se^2+Wh2se^2-2*Intercept*h2se*Wh2se)
  return(Z)
}

#
#
vars<-c("YearsEducation","BMI","Height","DrinksPerWeek","SevereObesity","AgeFirstBirth","BreastCancerFemale","PhysicalActivity","Type1Diabetes","HealthRating")
Label<-c("Years of Education","BMI","Height","Drinks per Week","Severe Obesity","Age at First Birth","Breast cancer",
         "Physical activity","Type 1 Diabetes","Self-rated health")

table<-matrix(ncol=13,nrow=length(vars))
table<-as.data.frame(table)
names(table)<-c("Phenotype","GWAS_h2","GWAS_h2_SE","IPWGWAS_h2","IPWGWAS_h2_SE","Intercept","P_dif","rG","rGSE","Intercepth2","Intercepth2SE","InterceptWh2","InterceptWh2SE")
table[,1]<-vars

count <- 1
for (pheno in vars) {
print(pheno)
path<-paste0("../OUTPUT/",pheno,"Results.h2.log") 
Wpath<-paste0("../OUTPUT/",pheno,"WeightedResults.h2.log") 
Rpath<-paste0("../OUTPUT/",pheno,".rg.log")

data<-readLines(path)
Wdata<-readLines(Wpath)
Rdata<-readLines(Rpath)

temp<-sub("h2:", "", data[grep("h2:",data)])   
h2<-as.numeric(strsplit(temp, " +")[[1]][4])
temp2<-strsplit(temp, " +")[[1]][5]
h2se<-as.numeric(substr(temp2,2,nchar(temp2)-1))
temp3<-sub("Intercept:", "", data[grep("Intercept:",data)]) 
Intercepth2<-as.numeric(strsplit(temp3, " +")[[1]][2])
temp4<-strsplit(temp3, " +")[[1]][3]
Intercepth2se<-as.numeric(substr(temp4,2,nchar(temp4)-1))

temp<-sub("h2:", "", Wdata[grep("h2:",Wdata)])   
Wh2<-as.numeric(strsplit(temp, " +")[[1]][4])
temp2<-strsplit(temp, " +")[[1]][5]
Wh2se<-as.numeric(substr(temp2,2,nchar(temp2)-1))
temp3<-sub("Intercept:", "", Wdata[grep("Intercept:",Wdata)]) 
InterceptWh2<-as.numeric(strsplit(temp3, " +")[[1]][2])
temp4<-strsplit(temp3, " +")[[1]][3]
InterceptWh2se<-as.numeric(substr(temp4,2,nchar(temp4)-1))

temp<-sub("Intercept:", "", Rdata[tail(grep("Intercept:",Rdata),n=1)]) 
Intercept<-as.numeric(strsplit(temp, " +")[[1]][2])

temp<-sub("Genetic Correlation:", "", Rdata[tail(grep("Genetic Correlation:",Rdata),n=1)]) 
rG<-as.numeric(strsplit(temp, "\\(")[[1]][1]) 
temp2<-strsplit(temp, "\\(")[[1]][2]
rGSE<-as.numeric(substr(temp2,2,nchar(temp2)-1)) 


Z<-PDifference(h2,Wh2,h2se,Wh2se,Intercept)
P<-2*pnorm(-abs(Z))

table[count,2]<-h2
table[count,3]<-h2se
table[count,4]<-Wh2
table[count,5]<-Wh2se
table[count,6]<-Intercept
table[count,7]<-P
table[count,8]<-rG
table[count,9]<-rGSE
table[count,10]<-Intercepth2
table[count,11]<-Intercepth2se
table[count,12]<-InterceptWh2
table[count,13]<-InterceptWh2se

count<-count+1
}

table<-cbind(Label,table)

table <- table[order(table$Label),]

head(table)

write.table(table,file="../TEMP/LDScoreResults.txt",sep="\t",quote=FALSE,row.names=FALSE)

table$GWAS_h2<-paste0(table$GWAS_h2," (",table$GWAS_h2_SE,")")
table$IPWGWAS_h2<-paste0(table$IPWGWAS_h2," (",table$IPWGWAS_h2_SE,")")
table$Intercepth2<-paste0(table$Intercepth2," (",table$Intercepth2SE,")")
table$InterceptWh2<-paste0(table$InterceptWh2," (",table$InterceptWh2SE,")")


pTemp<-table$P_dif
pTemp[table$P_dif<1e-09]<-paste0("$",substr(as.character(format(table$P_dif[table$P_dif<1e-09]),digits=16),1,4),"\\times 10^{-",substr(as.character(table$P_dif[table$P_dif<1e-09]),nchar(as.character(table$P_dif[table$P_dif<1e-09]))-1
                                                                                                                                                                               ,nchar(as.character(table$P_dif[table$P_dif<1e-09]))),"}$")
pTemp[table$P_dif>=1e-09&table$P_dif<=1e-03]<-paste0("$",substr(as.character(format(table$P_dif[table$P_dif>=1e-09&table$P_dif<=1e-03]),digits=16),1,4),"\\times 10^{-",substr(as.character(table$P_dif[table$P_dif>=1e-09&table$P_dif<=1e-03]),nchar(as.character(table$P_dif[table$P_dif>=1e-09&table$P_dif<=1e-03]))
                                                                                ,nchar(as.character(table$P_dif[table$P_dif>=1e-09&table$P_dif<=1e-03]))),"}$")

pTemp[table$P_dif>1e-03]<-round(table$P_dif[table$P_dif>1e-03],digits=3)

table$P_dif<-pTemp

table<-table[c("Label","GWAS_h2","IPWGWAS_h2","P_dif","Intercepth2","InterceptWh2")]

xTab<-xtable(table,align = "ll|l|l|l|l|l", caption="\\textbf{Caption here:} describe p-values")

colnames(xTab) = c("\\textbf{Phenotype}", "\\textbf{GWAS $h^2$} (SE)", "\\textbf{WGWAS $h^2$} (SE)","\\textbf{P}","\\textbf{GWAS Intercept} (SE)","\\textbf{WGWAS Intercept} (SE)" )
print(xTab,include.rownames=FALSE,sanitize.colnames.function = paste0,
      sanitize.text.function= paste0,type="latex",file="../OUTPUT/TABLES/LDScoreHeritabilities.tex",floating=FALSE)
