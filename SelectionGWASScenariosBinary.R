#Goal: Simulate some data to illustrate the effects of collider bias in the UKB:
library(tidyverse)
library(ggExtra)
library(stargazer)
library(xtable)

setwd("C:/Users/sjoerd/Dropbox/UKB_WGWAS/CODE")

set.seed(3454353)

for (p in c(0.05,0.25,0.5,0.75,0.95)){
N<-1000000
MAF<-0.4
SNP<-rbinom(N,2,MAF) #simulate SNP data with MAF of 0.4
table(SNP)
Error<-rnorm(N, mean = 0, sd = 1)

Y<-0.04*SNP+Error

#Make Y binary:
#p<-0.05

Y <- as.numeric(Y > quantile(Y,1-p))

#Scenario 1a SNP and y positively related, selection on y positive:
Data<-data.frame(SNP,Y,Error)

# Set the sampling probabilities based on the value of y
sampling_prob <- ifelse(Data$Y == 1, 2, 1)

# Normalize the sampling probabilities
sampling_prob <- sampling_prob / sum(sampling_prob)

# Sample from the original data with the specified sampling probabilities
data_oversampled <- Data[sample(1:0.05*N, size = 0.05*N, replace = TRUE, prob = sampling_prob), ]

sampling_prob <- ifelse(Data$Y == 1, 0.5, 1)
sampling_prob <- sampling_prob / sum(sampling_prob)
data_undersampled <- Data[sample(1:0.05*N, size = 0.05*N, replace = TRUE, prob = sampling_prob), ]


FillColors<-c("#E69F00", "#56B4E9")

LMTrue<-lm(Y~SNP,data=Data)
LMTrue$prev<-mean(Data$Y)
LMScen1_over<-lm(Y~SNP,data=data_oversampled)
LMScen1_over$prev<-mean(data_oversampled$Y)
LMScen1_under<-lm(Y~SNP,data=data_undersampled)
LMScen1_under$prev<-mean(data_undersampled$Y)

LMTrue$MAF<-mean(Data$SNP/2)

LMScen1_over$MAF<-mean(data_oversampled$SNP/2)
LMScen1_under$MAF<-mean(data_undersampled$SNP/2)

#Scenario 2a, SNP and y positivelY related, selection on SNP and y positive:
Data<-data.frame(SNP,Y,Error)

sampling_prob<-ifelse(Data$Y == 1 & Data$SNP == 2,p+0.04, 
               ifelse(Data$Y ==1 & Data$SNP == 1, p+0.03,
               ifelse(Data$Y ==0 & Data$SNP == 2, p+0.02,
               ifelse(Data$Y ==1 & Data$SNP == 0, p+0.02,
               ifelse(Data$Y ==0 & Data$SNP == 1, p+0.01,p)))))

sampling_prob <- sampling_prob / sum(sampling_prob)
data_oversampled <- Data[sample(1:0.05*N, size = 0.05*N, replace = TRUE, prob = sampling_prob), ]

sampling_prob<-ifelse(Data$Y == 1 & Data$SNP == 2,p-0.04, 
               ifelse(Data$Y ==1 & Data$SNP == 1, p-0.03,
               ifelse(Data$Y ==0 & Data$SNP == 2, p-0.02,
               ifelse(Data$Y ==1 & Data$SNP == 0, p-0.02,
               ifelse(Data$Y ==0 & Data$SNP == 1, p-0.01,p)))))

sampling_prob <- sampling_prob / sum(sampling_prob)

data_undersampled <- Data[sample(1:0.05*N, size = 0.05*N, replace = TRUE, prob = sampling_prob), ]

LMScen2a<-lm(Y~SNP,data=Data)
LMScen2a$prev<-mean(Data$Y)
LMScen2a_over<-lm(Y~SNP,data=data_oversampled)
LMScen2a_over$prev<-mean(data_oversampled$Y)

LMScen2a_under<-lm(Y~SNP,data=data_undersampled)
LMScen2a_under$prev<-mean(data_undersampled$Y)

LMScen2a_over$MAF<-mean(data_oversampled$SNP/2)
LMScen2a_under$MAF<-mean(data_undersampled$SNP/2)


#Scenario 2b, y positively, SNP negatively related, selection on SNP and y positive:
Data<-data.frame(SNP,Y,Error)

sampling_prob<-ifelse(Data$Y == 1 & Data$SNP == 2, p, 
                ifelse(Data$Y ==1 & Data$SNP == 1, p+0.01,
                ifelse(Data$Y ==0 & Data$SNP == 2, p-0.02,
                ifelse(Data$Y ==1 & Data$SNP == 0, p+0.02,
                ifelse(Data$Y ==0 & Data$SNP == 1, p-0.01,p)))))

sampling_prob <- (sampling_prob / sum(sampling_prob))*10000

data_oversampled <- Data[sample(1:0.05*N, size = 0.05*N, replace = TRUE, prob = sampling_prob), ]

sampling_prob<-ifelse(Data$Y == 1 & Data$SNP == 2, p, 
                      ifelse(Data$Y ==1 & Data$SNP == 1, p-0.01,
                             ifelse(Data$Y ==0 & Data$SNP == 2, p+0.02,
                                    ifelse(Data$Y ==1 & Data$SNP == 0, p-0.02,
                                           ifelse(Data$Y ==0 & Data$SNP == 1, p+0.01,p)))))

sampling_prob <- (sampling_prob / sum(sampling_prob))*10000


data_undersampled <- Data[sample(1:0.05*N, size = 0.05*N, replace = TRUE, prob = sampling_prob), ]

LMScen2b<-lm(Y~SNP,data=Data)
LMScen2b_over<-lm(Y~SNP,data=data_oversampled)
LMScen2b_over$prev<-mean(data_oversampled$Y)

LMScen2b_under<-lm(Y~SNP,data=data_undersampled)
LMScen2b_under$prev<-mean(data_undersampled$Y)


LMScen2b_over$MAF<-mean(data_oversampled$SNP/2)
LMScen2b_under$MAF<-mean(data_undersampled$SNP/2)


library(xtable)

# Function to extract necessary statistics from the model
extract_stats <- function(model) {
  n <- length(model$residuals)
  rsq <- sprintf("%.4f", summary(model)$r.squared)
  print(rsq)
  # Extract the coefficient and standard error for the SNP variable and the constant
  coef_snp <- round(coef(model)["SNP"],digits=4)
  coef_constant <- round(coef(model)["(Intercept)"],digits=4)
  se_snp <- round(summary(model)$coefficients["SNP", "Std. Error"],digits=5)
  se_constant <- round(summary(model)$coefficients["(Intercept)", "Std. Error"],digits=5)
  
  prev <- round(model[["prev"]],digits=5)
  maf <- round(model[["MAF"]],digits=2)
  return(c(paste0(coef_snp," (",se_snp,")"), paste0(coef_constant," (",se_constant,")"), n, rsq,prev,maf))
}

# Combine statistics from all models
stats_over <- cbind(
  extract_stats(LMTrue),
  extract_stats(LMScen1_over),
  extract_stats(LMScen2a_over),
  extract_stats(LMScen2b_over)
)

stats_under <- cbind(
  extract_stats(LMScen1_under),
  extract_stats(LMScen2a_under),
  extract_stats(LMScen2b_under)
)

# Name rows and columns for readability
colnames(stats_over) <- c("Population", "Scenario 1", "Scenario 2a", "Scenario 2b")
stats_under<-cbind(NA,stats_under)
colnames(stats_under) <- colnames(stats_over) 

rownames(stats_over) <- rownames(stats_under) <- c("SNP","Constant","Observations", "R-squared","prob([Y=1])","MAF")

# Convert stats matrices to data frames
stats_over_df <- as.data.frame(stats_over)
stats_under_df <- as.data.frame(stats_under)

stats_over_df$rownames<-rownames(stats_over_df)
stats_under_df$rownames<-rownames(stats_under_df)

combined_df <- rbind(stats_over_df, 
                     c(NA, NA, NA, NA, NA, NA), stats_under_df)

combined_df <- combined_df[, c(ncol(combined_df), 1:(ncol(combined_df) - 1))]
combined_df$Method <- ""
combined_df$Method[1] <- "Oversampling"
combined_df$Method[7] <- "Undersampling"

combined_df <- combined_df[, c(ncol(combined_df), 1:(ncol(combined_df) - 1))]

colnames(combined_df)[1] <- ""
colnames(combined_df)[2] <- ""
# Create a single LaTeX table using xtable
latex_table <- xtable(combined_df, caption = "Phenotype", align = c("l", "l","l", rep("c", ncol(combined_df) - 2)))

print(latex_table, type = "latex", append = FALSE, file = paste0("../OUTPUT/TABLES/SimulationResultsBinary",p,".tex"), 
      floating=FALSE,add.to.row=list(pos=list(-1,6),command = c("\\hline","\\hline ")), 
      include.rownames = FALSE)
}
