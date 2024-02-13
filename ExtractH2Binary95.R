# Load necessary libraries
library(ggplot2)
library(tidyr)

iter<-15

for (caus in c("1","10","200","2000","-1")){

h2mat1<-matrix(0, nrow = iter, ncol = 6)
h2mat2a<-matrix(0, nrow = iter, ncol = 6)
h2mat2b<-matrix(0, nrow = iter, ncol = 6)


j<-0
for (j in 1:6){
  for (i in 1:iter){

path<-paste0("../OUTPUT/Simulations/Simul1Binary95",i,j,"_",caus,".h2.log") 
data<-readLines(path)
temp<-sub("h2:", "", data[grep("h2:",data)])   
h2mat1[i,j]<-as.numeric(strsplit(temp, " +")[[1]][4])

path<-paste0("../OUTPUT/Simulations/Simul2aBinary95",i,j,"_",caus,".h2.log") 
data<-readLines(path)
temp<-sub("h2:", "", data[grep("h2:",data)])   
h2mat2a[i,j]<-as.numeric(strsplit(temp, " +")[[1]][4])

path<-paste0("../OUTPUT/Simulations/Simul2bBinary95",i,j,"_",caus,".h2.log") 
data<-readLines(path)
temp<-sub("h2:", "", data[grep("h2:",data)])   
h2mat2b[i,j]<-as.numeric(strsplit(temp, " +")[[1]][4])

  }
}


h2mat1<-apply(h2mat1,2,mean)
h2mat2a<-apply(h2mat2a,2,mean)
h2mat2b<-apply(h2mat2b,2,mean)
betalist<-c(0,0.2,0.4,0.6,0.8,1)

data<-cbind(h2mat1,h2mat2a,h2mat2b,betalist)

# Reshape the data into long format
data<-as.data.frame(data)
long_data <- pivot_longer(data, cols = -betalist, names_to = "variable", values_to = "value")
long_data$scenario <- factor(long_data$variable, 
                             levels = c("h2mat1", "h2mat2a", "h2mat2b"),
                             labels = c("scenario 1", "scenario 2a", "scenario 2b"))

# Plot the graph
plot<-ggplot(long_data, aes(x = betalist, y = value, color = scenario, linetype = scenario)) +
  geom_line() +
  geom_hline(yintercept = 0.2, linetype = "dashed") +
  scale_y_continuous(limits = c(min(long_data$value)-0.01, max(long_data$value)+0.01)) +
  labs(title = "",
       x = "Beta",
       y = "Heritability") +
  scale_x_continuous(breaks = c(0,0.2, 0.4, 0.6, 0.8, 1.0))+
  theme_minimal() +
  guides(color = guide_legend("Scenario"), linetype = guide_legend("Scenario"))

ggsave(plot=plot,height=6,width=6,dpi=200, filename=paste0("../OUTPUT/FIGURES/H2SimulationsBinary95",caus,".pdf"))

scen1<-read.table(paste0("../TEMP/Simul1ResultsBinary95",caus,".txt"),sep="\t")

scen1<-abs(scen1)
truth<-scen1$V7
scen1$V7<-NULL
#scen1_1<-scen1[1:6]
#scen1_2<-scen1[7:12]
#names(scen1_2)<-names(scen1_1)
#scen1<-rbind(scen1_1,scen1_2)

scen2a<-read.table(paste0("../TEMP/Simul1Results2aBinary95",caus,".txt"),sep="\t")
scen2b<-read.table(paste0("../TEMP/Simul1Results2bBinary95",caus,".txt"),sep="\t")

scenmat1<-apply(scen1,2,mean)
scenmat2a<-apply(scen2a,2,mean)
scenmat2b<-apply(scen2b,2,mean)


scenmat2a<-abs(scenmat2a)
scenmat2b<-abs(scenmat2b)

betalist<-c(0,0.2,0.4,0.6,0.8,1)

data<-cbind(scenmat1,scenmat2a,scenmat2b,betalist)


# Reshape the data into long format
data<-as.data.frame(data)
long_data <- pivot_longer(data, cols = -betalist, names_to = "variable", values_to = "value")
long_data$scenario <- factor(long_data$variable, 
                             levels = c("scenmat1", "scenmat2a", "scenmat2b"),
                             labels = c("scenario 1", "scenario 2a", "scenario 2b"))

# Plot the graph
plot<-ggplot(long_data, aes(x = betalist, y = value, color = scenario, linetype = scenario)) +
  geom_line() +
  geom_hline(yintercept = truth, linetype = "dashed") +
  scale_y_continuous(limits = c(min(long_data$value)-abs(max(long_data$value)), max(long_data$value)+abs(max(long_data$value)))) +
  labs(title = "",
       x = "Beta",
       y = "Point estimate") +
  scale_x_continuous(breaks = c(0,0.2, 0.4, 0.6, 0.8, 1.0))+
  theme_minimal() +
  guides(color = guide_legend("Scenario"), linetype = guide_legend("Scenario"))

ggsave(plot=plot,height=6,width=6,dpi=200, filename=paste0("../OUTPUT/FIGURES/PointEstSimulationsBinary95",caus,".pdf"))
}
