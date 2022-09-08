########################################################################################
## Description: Cancer detection on rare prostate morphologies (Test set 5) in CP MS

## Test data 1: Containing 27 slides scanned in Uppsala containing of special rare morphologies. 
# Scanning round #1 ('lars_special_cases')
# Including subtypes: Adenosis, Postatrophic hyperplasia, Atrophic-cancer, 
#                     PIN like cancer, Pseudohyper-plastic cancer

## Test data 2: Containing 152 slides scanned in Uppsala containing of special rare morphologies. 
# Scanning round #2 During summer 2022
# Including benign-subtypes: Adenosis, Basal_cell_hyperplasia, Clearcell_cribriform_hyperplasia
#                            Prostatic_atrophy, Postatrophic hyperplasia, Cowper
# Including cancer-subtypes: Atrophy_like_cancer, Foamy, Pseudohyperplastic, Small-cell-cancer

########################################################################################

rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("pROC")
source("functions.R")

########################################################################################
## Set Parameters
########################################################################################
use_testsplit <- FALSE
write_output <- FALSE 
datasource <- 'both'

########################################################################################
## Rare-subtypes
########################################################################################
# df1
df_test1 <- read.csv(file = file.path('Data', 'Demo', 
                                      'df_lars_special_cases_argmax.csv'))
# df2
df_test2 <- read.csv(file = file.path('Data', 'Demo',
                                      'rare_subtypes_ks.csv'))

########################################################################################
## Table Prediction-regions by subtype and Re-code subtypes
########################################################################################
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
## Combine Testset #1 & Testset #2
########################################################################################
keepvars <- c('pr_cx', 'pr_ben', "cl_slide_cx", 'ISUP', 'cx', 'subtype')
df_test1 <- select(df_test1, keepvars)
df_test2 <- select(df_test2, keepvars)

########################################################################################
## Stratify on different data sources
########################################################################################
if (datasource == 'test1'){
  df_test <- df_test1
} else if (datasource == 'test2'){
  df_test <- df_test2
} else {
  df_test <- rbind(df_test1, df_test2)
}

########################################################################################
## Sample part of Internal Test-Data as calibration
########################################################################################
if (use_testsplit){
  # II. Add test-split from I.
  df_calib <- read.csv(file = file.path('Data', 'Demo', 'df_test.csv'))
  df_calib <- filter(df_calib, calibset == 1)
} else {
  df_calib <- read.csv(file = file.path('Data', 'Demo', 'df_test.csv'))
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

########################################################################################
## Reorder columns and save output
########################################################################################
reorder_columns <- function(tab){
  x <- tab
  df <- as.data.frame.matrix(tab) 
  columnorder <- c('Adenosis', 'Basal_cell_hyperplasia', 'Clearcell_cribriform_hyperplasia',
                   'Prostatic_atrophy', 'Postatrophic_hyperplasia', 'Cowper',
                   'Atrophy_like_cancer', 'Foamy', 'PIN_like_cx', 
                   'Pseudohyperplastic', 'Small.cell.cancer', 'Sum')
  
  mytable <- df[, columnorder]
  
  return(mytable)
}
tab01 <- reorder_columns(tab01)
tab05 <- reorder_columns(tab05)
tab1 <- reorder_columns(tab1)

if (write_output){
  write.table(tab01, file = file.path('Output', 'Tables', 'CP_CX_rare_subtype.csv'))
  write.table(tab05, file = file.path('Output', 'Tables', 'CP_CX_rare_subtype.csv'), append = TRUE)
  write.table(tab1, file = file.path('Output', 'Tables', 'CP_CX_rare_subtype.csv'), append = TRUE)
}
tab01;tab05;tab1

################################## end of program #################################
