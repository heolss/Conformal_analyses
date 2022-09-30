
rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
source("functions.R")

## Original Data
df_test <- read.csv(file = file.path('Data', 'Demo', 'df_test.csv'))

########################################################################################
## Set Parameters
########################################################################################
normalize <- FALSE
create_testsplit <- FALSE
write_output <- FALSE 

set.seed(1)

########################################################################################
## Sample part of the training data as calibration set
# Create split based on man-level (unique individuals and all their corresponding biopsies
########################################################################################
# I. Run only once
if (create_testsplit){
  ids <- sample(unique(df_test$man), 123)
  df_calib <- df_test[df_test$man %in% ids, ]
  df_test$calibset <- ifelse(df_test$man %in% ids, 0, 1)
  write.csv(df_test, file = file.path('Data', 'Demo', 'df_test.csv'))
}

# II. Add test-split from I.
df_calib <- filter(df_test, calibset == 1)
df_test <- filter(df_test, calibset == 0)

########################################################################################
## Create CV-folds
########################################################################################
## Original Data
a <- df_test

# randomly shuffle rows
rows <- sample(nrow(df_test))
df_test <- df_test[rows, ]

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
df <- filter(df, testLabels == 2)
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

BASECASE_detect_drift <- ggplot(data=df2, aes(x=len, y=p_ks)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0, 100, 10)) +   
  ylab("P-value") +
  xlab("Percent of data") + 
  # ggtitle("Testset 1") + 
  geom_hline(yintercept = .05, linetype="dashed", color = "red")


if (write_output){
  ggsave(filename = file.path('Output', 'Figures', 'BASECASE_detect_drift.png'), 
         BASECASE_detect_drift)
}

################################## end of program #################################
