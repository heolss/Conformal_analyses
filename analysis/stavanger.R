########################################################################################
## Description:  
# Analysis - External scanner and external pathology laboratory - Test set 5
# Author Henrik Olsson
########################################################################################

rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
source("functions.R")

## Original Data
df_test <- read.csv(file = file.path('Data', 'Demo', 'stavanger.csv'))

########################################################################################
## Set Parameters
########################################################################################
normalize <- FALSE
write_output <- FALSE 
use_testsplit <- FALSE

########################################################################################
## Sample part of Internal Test-Data as calibration

# Either use the full Test-data as calibrationset or the same test-split created for CP_LO
########################################################################################
if (use_testsplit){
  # II. Add test-split from I.
  df_calib <- read.csv(file = file.path('Data', 'Demo', 'df_test.csv'))
  df_calib <- filter(df_calib, calibset == 1)
} else {
  df_calib <- read.csv(file = file.path('Data', 'Demo', 'df_test.csv'))
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
cx <- calib_plot(pValues_cx, testSet = testSet)
# Cross-tabulate Outcome CX V. CP-CX Or Outcome ISUP V. CP-ISUP
df_pred01 <- df_predict_region(pValues_cx, sigfLevel = .001, outcome = 'CX')
df_pred05 <- df_predict_region(pValues_cx, sigfLevel = .005, outcome = 'CX')
df_pred1 <- df_predict_region(pValues_cx, sigfLevel = .01, outcome = 'CX')
tab01 <- tab_predict_region(df_pred01)
tab05 <- tab_predict_region(df_pred05)
tab1 <- tab_predict_region(df_pred1)
if (write_output){
  write.table(tab01, file = file.path('Output', 'Tables', 'stavanger_CX.csv'))
  write.table(tab05, file = file.path('Output', 'Tables', 'stavanger_CX.csv'), append = TRUE)
  write.table(tab1, file = file.path('Output', 'Tables', 'stavanger_CX.csv'), append = TRUE)
  # Cancer-detection plot
  cx <- cx + theme(text = element_text(size = 17))
  cx <- cx + ggtitle("")
  ggsave(filename = file.path('Output', 'Figures', 'cx_stavanger.png'), cx, 
         width = 6, 
         height = 6)
}
tab01;tab05;tab1


################################## end of program #################################
