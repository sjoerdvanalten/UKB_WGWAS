library(ggplot2)

vars<-c("EA3","CigarettesPerDay","DrinksPerWeek","SmokingInitiation","BMI2",
        "HIP","WHR","WC","Height","Schizophrenia",
        "Depression","DBP","SBP","HTN","Allergic",
        "SWB","AFB","ParticipationFoodQuestionnaire","ParticipationPhysicalMonitor","ParticipationMentalHealthQuestionnaire",
        "Type1Diabetes","BreastCancer")

labels<-c("Educational Attainment","Cigarettes per Day","Drinks per Week","Smoking Initiation","BMI",
          "Hip Circumference","Waist to Hip Ratio","Waist Circumference","Height","Schizophrenia",
          "Depression","Diastollic blood pressure","Systollic blood pressure","Hypertension","Allergies",
          "Subjective Well-being","Age at First Birth","Participation (food questionnaire, UKB)","Participation (physical monitor, UKB)",
          "Participation (mental health questionnaire, UKB)","Diabetes (type 1)","Breast Cancer")

table<-matrix(ncol=4,nrow=length(vars))
table<-as.data.frame(table)
names(table)<-c("Phenotype","rG","se","p")
table[,1]<-vars

count<-1

for (pheno in vars) {
path<-paste0("../OUTPUT/GWASWeightsand",pheno,".rg.log")
data<-readLines(path)
index<-grep("p1",data)+1
tmp<-strsplit(data[index], " ")
tmp<-tmp[[1]]
tmp<-tmp[tmp!=""]
table$rG[count]<-as.numeric(tmp[3])
table$se[count]<-as.numeric(tmp[4])
table$p[count]<-as.numeric(tmp[6])
count<-count+1
}

alpha <- 0.05
bonferroni_alpha <- alpha / 22
z_score <- qnorm(1 - bonferroni_alpha / 2)

table$lower<- table$rG - z_score * table$se
table$higher <- table$rG + z_score * table$se



table<-cbind(table,labels)

ggplot(table,aes(x=reorder(labels,-rG),y=rG,ymin=lower,ymax=higher))+
  geom_point()+
  geom_errorbar()+
  coord_flip()+
  theme_bw() + 
  scale_y_continuous(limits=c(-0.8,0.8)) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") + xlab("rG")
ggsave("../OUTPUT/RGGWASOnTheWeights.pdf")