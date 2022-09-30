rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
source("functions.R")

## Original Data
df_test <- read.csv(file = file.path('Data', 'Demo', 'sthlm3_aperio_hamamatsu.csv'))

########################################################################################
## Set Parameters
########################################################################################
write_output <- FALSE 
set.seed(1)

########################################################################################
## Add STHLM3_APERIO Calibration set
########################################################################################
df_calib <- read.csv(file = file.path('Data', 'Demo', 'sthlm3_aperio_calib.csv'))

########################################################################################
## Only use Aperio slides of Test-Set as CalibrationSet
########################################################################################
df_calib <- filter(df_calib, ext == '.svs')
table(df_calib$ext)

########################################################################################
## Create CV-folds
########################################################################################
a <- df_test

#Create k equally size folds
folds <- cut(seq(1,nrow(a)),breaks=100,labels=FALSE)

########################################################################################
## Calculate p-value on first fold
########################################################################################
df_total <- data.frame()
p_ks <- data.frame()

testIndexes <- which(folds==1,arr.ind=TRUE)
df_test <- a[testIndexes, ]

## Generate prediction-regions for CX
data <- generate_vars(outcome = "CX", testdata = "test")
testSet <- data$testSet
testLabels <- data$testLabels
calibset <- data$df_calib
col_outcome <- data$col_outcome
# ICP P-values: Class-conditional Inductive conformal classifier for multi-class problems
pValues_cx <- ICPClassification(testSet = testSet, calibSet = calibset)
## P-value
df <- data.frame(testLabels, pValues_cx)
# df <- filter(df, testLabels == 2)
# Only unique - ties not allowed
df <- df[!duplicated(df$X2),]
p_ks <- ks.test(df$X2, 'punif', 0, 1)$p.value

########################################################################################
## Loop over remaining folds
########################################################################################
for (i in 2:100){
  # Create 
  testIndexes <- which(folds %in% c(i),arr.ind=TRUE)
  df <- a[testIndexes, ]
  df_test <- rbind(df_test, df)
  
  ## Generate prediction-regions for CX
  data <- generate_vars(outcome = "CX", testdata = "test")
  testSet <- data$testSet
  testLabels <- data$testLabels
  calibset <- data$df_calib
  col_outcome <- data$col_outcome
  # ICP P-values: Class-conditional Inductive conformal classifier for multi-class problems
  pValues_cx <- ICPClassification(testSet = testSet, calibSet = calibset)
  ## P-value
  df <- data.frame(testLabels, pValues_cx)
  # Only unique - ties not allowed
  df <- df[!duplicated(df$X2),]
  p <- ks.test(df$X2, 'punif', 0, 1)$p.value
  p_ks <- rbind(p_ks, p)
}
p_ks

########################################################################################
## Lineplot of results
# Kolmogorov P-value of equality in distribution of test data with training data
########################################################################################
len <- 1:100
df2 <- data.frame(p_ks, len)

HAMA_detect_drift <- ggplot(data=df2, aes(x=len, y=p_ks)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 100, 5)) +   
  ylab("P-value") +
  xlab("Percent of data") + 
  geom_hline(yintercept = .05, linetype="dashed", color = "red")

if (write_output){
  ggsave(filename = file.path('Output', 'Figures', 'HAMA_detect_drift.png'), 
         HAMA_detect_drift)
}

################################## end of program #################################
