
rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
source("Programs/functions.R")

## Original Data
df_test <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_imagebase.csv'))
df_ks1 <- read.csv(file = file.path('Data', 'Derived', 'conformal_20200225', 'df_KS1.csv'))
df_KS1_CGAN <- read.csv(file = file.path('Data', 'Derived', 'conformal_20200225', 'df_KS1_CGAN.csv'))

# Evaluate against the mode ISUP assessment of the 23 pathologists in the Imagebase-panel
imagebase <- read.csv(file = file.path('Data', 'Derived', 'CP_LO', 'FinalResultsProstate_with_AI_anonymized.csv'))
df_test <- left_join(df_test, select(imagebase, c(Path.no, MODE, Consensus)), by = c("slide" = "Path.no"))
df_test <- filter(df_test, !is.na(MODE))
df_test$ISUP <- df_test$MODE
df_test <- df_test[order(df_test$slide),]

########################################################################################
## Set Parameters
########################################################################################
normalize <- FALSE
ratioTrain = 0.5
write_output <- FALSE 
use_testsplit <- FALSE

########################################################################################
## Sample part of the training data as calibration set
########################################################################################
if (use_testsplit){
  # II. Add test-split from I.
  df_calib <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_test.csv'))
  df_calib <- filter(df_calib, calibset == 1)
} else {
  df_calib <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_test.csv'))
}

########################################################################################
## Normalize predictions:
########################################################################################
if (normalize){
  cols <- c('X0','X1','X2','X3','X4','X5')
  r_sum <- rowSums(df_calib[cols])
  df_calib$X0 <- df_calib$X0 / r_sum
  df_calib$X1 <- df_calib$X1 / r_sum
  df_calib$X2 <- df_calib$X2 / r_sum
  df_calib$X3 <- df_calib$X3 / r_sum
  df_calib$X4 <- df_calib$X4 / r_sum
  df_calib$X5 <- df_calib$X5 / r_sum
  cols <- c('X0','X1','X2','X3','X4','X5')
  r_sum <- rowSums(df_test[cols])
  df_test$X0 <- df_test$X0 / r_sum
  df_test$X1 <- df_test$X1 / r_sum
  df_test$X2 <- df_test$X2 / r_sum
  df_test$X3 <- df_test$X3 / r_sum
  df_test$X4 <- df_test$X4 / r_sum
  df_test$X5 <- df_test$X5 / r_sum
}

########################################################################################
## Confusion matrix and Accuracy
########################################################################################
caret::confusionMatrix(factor(df_calib$cl_slide_isup), factor(df_calib$ISUP)); vcd::Kappa(table(df_calib$ISUP, df_calib$cl_slide_isup))
caret::confusionMatrix(factor(df_test$cl_slide_isup), factor(df_test$ISUP)); vcd::Kappa(table(df_test$ISUP, df_test$cl_slide_isup))

########################################################################################
## Generate Output for TestData
########################################################################################
## CX
# Generate Variables: Test-Data / Calibration-Set / Outcome
data <- generate_vars(outcome = "CX", testdata = "test")
testSet <- data$testSet
testLabels <- data$testLabels
calibset <- data$df_calib
col_outcome <- data$col_outcome
# ICP P-values: Class-conditional Inductive conformal classifier for multi-class problems
pValues_cx <- ICPClassification(testSet = testSet, calibSet = calibset)
# Calibration Plot
cx <- calib_plot(pValues_cx, testSet = testSet) + ggtitle("Imagebase: CX")
# Compute conformal prediction regions at different significance levels
df_pred01 <- df_predict_region(pValues_cx, sigfLevel = .001, outcome = 'CX')
df_pred05 <- df_predict_region(pValues_cx, sigfLevel = .005, outcome = 'CX')
df_pred1 <- df_predict_region(pValues_cx, sigfLevel = .01, outcome = 'CX')
# Cross-tabulate Outcome CX V. CP-CX Or Outcome ISUP V. CP-ISUP
tab01 <- tab_predict_region(df_pred01)
tab05 <- tab_predict_region(df_pred05)
tab1 <- tab_predict_region(df_pred1)
if (write_output){
  write.table(tab01, file = file.path('Output', 'Tables', 'PredictionSets', 'iamgebase_CX.csv'))
  write.table(tab05, file = file.path('Output', 'Tables', 'PredictionSets', 'iamgebase_CX.csv'), append = TRUE)
  write.table(tab1, file = file.path('Output', 'Tables', 'PredictionSets', 'iamgebase_CX.csv'), append = TRUE)
}
tab01;tab05;tab1

## ISUP
# Generate Variables: Test-Data / Calibration-Set / Outcome
data <- generate_vars(outcome = "ISUP", testdata = "test")
testSet <- data$testSet
testLabels <- data$testLabels
calibset <- data$df_calib
col_outcome <- data$col_outcome
# ICP P-values: Class-conditional Inductive conformal classifier for multi-class problems
pValues_isup <- ICPClassification(testSet = testSet, calibSet = calibset)
# Calibration Plot
isup <- calib_plot(pValues_isup, testSet = testSet) + ggtitle("Imagebase: ISUP")
# Compute conformal prediction regions at different significance levels
df_pred05 <- df_predict_region(pValues_isup, sigfLevel = .05, outcome = 'ISUP')
df_pred10 <- df_predict_region(pValues_isup, sigfLevel = .1, outcome = 'ISUP')
df_pred20 <- df_predict_region(pValues_isup, sigfLevel = .2, outcome = 'ISUP')
df_pred33 <- df_predict_region(pValues_isup, sigfLevel = 1/3, outcome = 'ISUP')
# Cross-tabulate Outcome CX V. CP-CX Or Outcome ISUP V. CP-ISUP
tab05 <- tab_predict_region(df_pred05)
tab10 <- tab_predict_region(df_pred10)
tab20 <- tab_predict_region(df_pred20)
tab33 <- tab_predict_region(df_pred33)
if (write_output){
  write.table(tab05, file = file.path('Output', 'Tables', 'PredictionSets', 'iamgebase_ISUP.csv'))
  write.table(tab10, file = file.path('Output', 'Tables', 'PredictionSets', 'iamgebase_ISUP.csv'), append = TRUE)
  write.table(tab20, file = file.path('Output', 'Tables', 'PredictionSets', 'iamgebase_ISUP.csv'), append = TRUE)
  write.table(tab33, file = file.path('Output', 'Tables', 'PredictionSets', 'iamgebase_ISUP.csv'), append = TRUE)
}
tab05;tab10;tab20;tab33

grid.arrange(cx, isup, ncol = 2)

########################################################################################
# Generate Table of predgroups by ISUP (Among CP-Multi-set-predictions
########################################################################################
multi10 <- predgroups(df_pred10)
multi20 <- predgroups(df_pred20)
multi33 <- predgroups(df_pred33)

if (write_output){
  write.table(multi10, file = file.path('Output', 'Tables', 'PredictionSets', 'multiset_IMB.csv'))
  write.table(multi20, file = file.path('Output', 'Tables', 'PredictionSets', 'multiset_IMB.csv'), append = TRUE)
  write.table(multi33, file = file.path('Output', 'Tables', 'PredictionSets', 'multiset_IMB.csv'), append = TRUE)
}

################################## end of program #################################
