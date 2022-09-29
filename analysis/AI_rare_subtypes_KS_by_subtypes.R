########################################################################################
## Description: Cancer detection on rare prostate morphologies (Test set 6) in CP MS

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
## Data-mgmt Testset #1
# Rename columns to match code and functions
# Add ISUP and CX from clinical database (database_20190130.csv)
########################################################################################
df_test1 <- read.csv(file = file.path('Data', 'Demo', 
                                      'df_lars_special_cases_argmax.csv'))

#######################################################################################
## Data-mgmt: Testset #2
########################################################################################
# AI-predictions
df_test2 <- read.csv(file = file.path('Data', 'Demo',
                                      'rare_subtypes_ks.csv'))

########################################################################################
## Table Prediction-regions by subtype and Re-code subtypes
########################################################################################
# Table helper function
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
keepvars <- c('pr_cx', 'pr_ben', 'ISUP', 'cx', 'subtype', 'cl_slide_cx')
df_test1 <- df_test1[keepvars]
df_test2 <- df_test2[keepvars]

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
## Table Point predictions (classification) by AI without CP
########################################################################################
# True/False classification by AI Lars-subtypes
df_test$correct_cx <- NULL
df_test <- mutate(df_test, correct_cx = ifelse(cx == cl_slide_cx, TRUE, FALSE))
table(df_test$correct_cx)

# Table
table_subtype <- function(df){
  x <- df
  tab <- table(x$correct_cx, x$subtype, useNA="ifany")
  tab_margins <- addmargins(tab, 2)
  mytable <- addmargins(tab)
  for (i in 1: dim(mytable)[1]-1){
    mytable[i,] <- paste(mytable[i,], " (", round( prop.table(tab_margins,2)[i,],2), ")", sep="")
  }
  return(mytable)
}

AI_CX_subtype <- table_subtype(df_test)

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
AI_CX_subtype <- reorder_columns(AI_CX_subtype)
AI_CX_subtype

########################################################################################
## Generate Output for TestData
########################################################################################
if (write_output){
  write.table(AI_CX_subtype, file = file.path('Output', 'Tables', 'AI_CX_rare_subtype.csv'))
}

################################## end of program #################################


