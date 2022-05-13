########################################################################################
## Description: lars_special_cases contains a dataset of 28 slides scanned in Uppsala
#               containing special rare morphologies. 

#               Correct and error predictions by the AI standalone without CP
########################################################################################

rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
source("Programs/functions.R")

## Original Data
# df_calib <- read.csv(file = file.path('Data', 'Derived', 'CP_LO', 'df_test.csv'))
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
AI_CX_subtype

########################################################################################
## Generate Output for TestData
########################################################################################
write.table(AI_CX_subtype, file = file.path('Output', 'Tables', 'PredictionSets', 'AI_CX_subtype.csv'))

################################## end of program #################################

