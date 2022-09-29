########################################################################################
## Description:  Analysis and Results - Baseline test set - Test-set 1
# Author Henrik Olsson
########################################################################################

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
## Generate Multi-label Output for TestData
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
cx <- calib_plot(pValues_cx, testSet = testSet) + ggtitle("Test: CX")
# Compute conformal prediction regions at different significance levels
df_pred01 <- df_predict_region(pValues_cx, sigfLevel = .001, outcome = 'CX')
df_pred05 <- df_predict_region(pValues_cx, sigfLevel = .005, outcome = 'CX')
df_pred1 <- df_predict_region(pValues_cx, sigfLevel = .01, outcome = 'CX')
# Cross-tabulate Outcome CX V. CP-CX Or Outcome ISUP V. CP-ISUP
tab01 <- tab_predict_region(df_pred01)
tab05 <- tab_predict_region(df_pred05)
tab1 <- tab_predict_region(df_pred1)
if (write_output){
  write.table(tab01, file = file.path('Output', 'Tables', 'test_CX.csv'))
  write.table(tab05, file = file.path('Output', 'Tables', 'test_CX.csv'), append = TRUE)
  write.table(tab1, file = file.path('Output', 'Tables', 'test_CX.csv'), append = TRUE)
  # Cancer-detection plot
  cx <- cx + theme(text = element_text(size = 17))
  cx <- cx + ggtitle("")
  ggsave(filename = file.path('Output', 'Figures', 'cx_test.png'), cx, 
         width = 6, 
         height = 6)
}
tab01;tab05;tab1

## ISUP
# Only ISUP1-ISUP5
df_test <- filter(df_test, ISUP != 0)
df_calib <- filter(df_calib, ISUP != 0)
# Generate Variables: Test-Data / Calibration-Set / Outcome
data <- generate_vars(outcome = "ISUP", testdata = "test")
testSet <- data$testSet
testLabels <- data$testLabels
calibset <- data$df_calib
col_outcome <- data$col_outcome
# ICP P-values: Class-conditional Inductive conformal classifier for multi-class problems
pValues_isup <- ICPClassification(testSet = testSet, calibSet = calibset)
# Calibration Plot
isup <- calib_plot(pValues_isup, testSet = testSet) + ggtitle("Test: ISUP")
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
  write.table(tab05, file = file.path('Output', 'Tables', 'test_ISUP.csv'))
  write.table(tab10, file = file.path('Output', 'Tables', 'test_ISUP.csv'), append = TRUE)
  write.table(tab20, file = file.path('Output', 'Tables', 'test_ISUP.csv'), append = TRUE)
  write.table(tab33, file = file.path('Output', 'Tables', 'test_ISUP.csv'), append = TRUE)
}
tab05;tab10;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

########################################################################################
# Generate Table of predgroups by ISUP (Among CP-Multi-set-predictions
########################################################################################
predgroups <- function(df_pred){
  x <- df_pred
  # x <- filter(x, pred_group == "Multiple") # OBS! if only filtering out multiple predictions
  x$true_cols <- gsub("[^0-9 ]", "", x$true_cols)
  x$true_cols <- str_replace(x$true_cols, "2", "ISUP1")
  x$true_cols <- str_replace(x$true_cols, "3", "ISUP2")
  x$true_cols <- str_replace(x$true_cols, "4", "ISUP3")
  x$true_cols <- str_replace(x$true_cols, "5", "ISUP4")
  x$true_cols <- str_replace(x$true_cols, "6", "ISUP5")
  
  # Tabulate prediction-regions by ISUP
  tab <- table(x$true_cols, x$testLabels, useNA="ifany")
  tab_margins <- addmargins(tab, 2)
  mytable <- addmargins(tab)
  for (i in 1: dim(mytable)[1]-1){
    mytable[i,] <- paste(mytable[i,], " (", round( prop.table(tab_margins,2)[i,],2), ")", sep="")
  }
  return(mytable)
}

predgroup20 <- predgroups(df_pred20)
predgroup33 <- predgroups(df_pred33)

if (write_output){
  write.table(predgroup20, file = file.path('Output', 'Tables', 'predgroups_TEST.csv'))
  write.table(predgroup33, file = file.path('Output', 'Tables', 'predgroups_TEST.csv'), append = TRUE)
}

################################## end of program #################################
