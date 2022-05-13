########################################################################################
## Description: Special rare morphology cases from the STHLM3 study pathology assessment 
#               containing a set of 27 slides scanned in Uppsala
########################################################################################

rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
source("Programs/functions.R")

## Original Data
df_test <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_lars_special_cases_argmax.csv'))
df_ks1 <- read.csv(file = file.path('Data', 'Derived', 'conformal_20200225', 'df_KS1.csv'))
df_KS1_CGAN <- read.csv(file = file.path('Data', 'Derived', 'conformal_20200225', 'df_KS1_CGAN.csv'))

########################################################################################
## Set Parameters
########################################################################################
normalize <- FALSE
ratioTrain = 0.5
write_output <- FALSE 
use_testsplit <- FALSE

########################################################################################
## Recode error in clinical data
# OBS! These needs to be checked with Lars
########################################################################################
table(df_test$cx); table(df_test$ISUP)
df_test$ISUP[df_test$ISUP == 5] <- 1
df_test$ISUP[df_test$slide == 'SP845_2012 1'] <- 1
df_test$ISUP[df_test$slide == 'SP1475_2011 6'] <- 1
df_test$ISUP[df_test$slide == 'SP1034_2012 1V'] <- 1
# df_test$ISUP[df_test$slide == 'SP960_2013 1A'] <- 1
df_test$comment <- as.character(df_test$comment)
df_test$comment[df_test$slide == 'SP960_2013 1A'] <- "PAH"
df_test$cx <- ifelse(df_test$ISUP >= 1, 1, 0)
table(df_test$cx); table(df_test$ISUP)

# Remove RP 
df_test <- filter(df_test, slide != 'SP1034_2012 1V')

########################################################################################
## Sample part of Internal Test-Data as calibration
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
## Table Prediction-regions by subtype
########################################################################################
# coding subtypes
df_test <- within(df_test, {
  subtype <- NA
  subtype[comment == 'PAH'] <- 1
  subtype[comment == 'adenosis'] <- 2
  subtype[comment == 'cancer of atrophic type'] <- 3
  subtype[comment == 'PIN like cx'] <- 4
  subtype[comment == 'pseudohyperplastic prostatic cancer'] <- 5
  subtype <- factor(subtype, levels = c(1,2,3,4,5), 
                    labels = c('Postatrophic hyperplasia', "Adenosis", "Atrophic_CX", "PIN_CX", "PH_CX"))
})

# Table by subtype
table_subtype <- function(df){
  x <- df
  tab <- table(x$pred_group, x$subtype, useNA="ifany")
  tab_margins <- addmargins(tab, 2)
  mytable <- addmargins(tab)
  for (i in 1: dim(mytable)[1]-1){
    mytable[i,] <- paste(mytable[i,], " (", round( prop.table(tab_margins,2)[i,],2), ")", sep="")
  }
  return(mytable)
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
cx <- calib_plot(pValues_cx, testSet = testSet) + ggtitle("Imagebase: CX")
# Cross-tabulate Outcome CX V. CP-CX Or Outcome ISUP V. CP-ISUP
df_pred01 <- df_predict_region(pValues_cx, sigfLevel = .001, outcome = 'CX')
df_pred05 <- df_predict_region(pValues_cx, sigfLevel = .005, outcome = 'CX')
df_pred1 <- df_predict_region(pValues_cx, sigfLevel = .01, outcome = 'CX')
tab01 <- table_subtype(cbind(df_pred01, select(df_test, subtype)))
tab05 <- table_subtype(cbind(df_pred05, select(df_test, subtype)))
tab1 <- table_subtype(cbind(df_pred1, select(df_test, subtype)))

if (write_output){
  write.table(tab01, file = file.path('Output', 'Tables', 'PredictionSets', 'lars_CX_subtype.csv'))
  write.table(tab05, file = file.path('Output', 'Tables', 'PredictionSets', 'lars_CX_subtype.csv'), append = TRUE)
  write.table(tab1, file = file.path('Output', 'Tables', 'PredictionSets', 'lars_CX_subtype.csv'), append = TRUE)
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
isup <- calib_plot(pValues_isup, testSet = testSet) + ggtitle("Imagebase: ISUP")
# Cross-tabulate Outcome CX V. CP-CX Or Outcome ISUP V. CP-ISUP
df_pred10 <- df_predict_region(pValues_isup, sigfLevel = .1, outcome = 'ISUP')
df_pred20 <- df_predict_region(pValues_isup, sigfLevel = .2, outcome = 'ISUP')
df_pred33 <- df_predict_region(pValues_isup, sigfLevel = 1/3, outcome = 'ISUP')
tab10 <- table_subtype(cbind(df_pred10, select(df_test, subtype)))
tab20 <- table_subtype(cbind(df_pred20, select(df_test, subtype)))
tab33 <- table_subtype(cbind(df_pred33, select(df_test, subtype)))
if (write_output){
  write.table(tab10, file = file.path('Output', 'Tables', 'PredictionSets', 'lars_ISUP_subtype.csv'))
  write.table(tab20, file = file.path('Output', 'Tables', 'PredictionSets', 'lars_ISUP_subtype.csv'), append = TRUE)
  write.table(tab33, file = file.path('Output', 'Tables', 'PredictionSets', 'lars_ISUP_subtype.csv'), append = TRUE)
}
tab10;tab20;tab33

grid.arrange(cx, isup, ncol = 2)


################################## end of program #################################
