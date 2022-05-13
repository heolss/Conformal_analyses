rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
source("Programs/functions.R")

## Original Data
# df_calib <- read.csv(file = file.path('/Users/henrikolsson1/Desktop/Sandbox/Transport/CP-LO', 'df_train.csv'))
df_ks1 <- read.csv(file = file.path('Data', 'Derived', 'conformal_20200225', 'df_KS1.csv'))
df_KS1_CGAN <- read.csv(file = file.path('Data', 'Derived', 'conformal_20200225', 'df_KS1_CGAN.csv'))
df_test <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_test.csv'))

########################################################################################
## Test
########################################################################################
# II. Add test-split from I.
df_test <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_test.csv'))
df_calib <- filter(df_test, calibset == 1)
df_test <- filter(df_test, calibset == 0)

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
p_test <- ks.test(df$X2, 'punif', 0, 1)$p.value

########################################################################################
## Imagebase
########################################################################################
## Original Data
df_test <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_imagebase.csv'))
df_ks1 <- read.csv(file = file.path('Data', 'Derived', 'conformal_20200225', 'df_KS1.csv'))
df_KS1_CGAN <- read.csv(file = file.path('Data', 'Derived', 'conformal_20200225', 'df_KS1_CGAN.csv'))
# Mode instead of LE
imagebase <- read.csv(file = file.path('Data', 'Derived', 'CP_LO', 'FinalResultsProstate_with_AI_anonymized.csv'))
df_test <- left_join(df_test, select(imagebase, c(Path.no, MODE, Consensus)), by = c("slide" = "Path.no"))
df_test <- filter(df_test, !is.na(MODE))
df_test$ISUP <- df_test$MODE
df_test <- df_test[order(df_test$slide),]
# Calibration-set
df_calib <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_test.csv'))

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
p_imagebase <- ks.test(df$X2, 'punif', 0, 1)$p.value

########################################################################################
## KS-Data
########################################################################################
# II. Add test-split from I.
df_test <- read.csv(file = file.path('Data', 'Derived', 'CP_LO', 'df_ks.csv'))

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
## Lars-special-cases: 28 slides scanned in Uppsala containing special rare morphologies
########################################################################################
df_test <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_lars_special_cases_argmax.csv'))
## Recode error in clinical data
# OBS! These needs to be checked with Lars
df_test$ISUP[df_test$ISUP == 5] <- 1
df_test$ISUP[df_test$slide == 'SP845_2012 1'] <- 1
df_test$ISUP[df_test$slide == 'SP1475_2011 6'] <- 1
df_test$ISUP[df_test$slide == 'SP1034_2012 1V'] <- 1
# df_test$ISUP[df_test$slide == 'SP960_2013 1A'] <- 1
df_test$comment <- as.character(df_test$comment)
df_test$comment[df_test$slide == 'SP960_2013 1A'] <- "PAH"
df_test$cx <- ifelse(df_test$ISUP >= 1, 1, 0)
# Calibration-set
df_calib <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_test.csv'))

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
p_lars <- ks.test(df$X2, 'punif', 0, 1)$p.value

########################################################################################
## STHLM3-Paired: STHLM3-Aperio: Paired slides scanned on Aperio
# Separate predictions on Calibrationset
########################################################################################
df_calib <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'df_test.csv'))
df_test <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'predictions_slide_rescanned_aperio.csv'))

## Data-mgmt
df_test <- rename(df_test, X0=slide_pred_isup_0, X1=slide_pred_isup_1,
                  X2=slide_pred_isup_2, X3=slide_pred_isup_3,
                  X4=slide_pred_isup_4, X5=slide_pred_isup_5,
                  pr_cx = slide_pred_cancer)
df_test$pr_ben = 1 - df_test$pr_cx
# Rename columns Calib-Set
df_calib <- rename(df_calib, X0=slide_pred_isup_0, X1=slide_pred_isup_1,
                   X2=slide_pred_isup_2, X3=slide_pred_isup_3,
                   X4=slide_pred_isup_4, X5=slide_pred_isup_5,
                   pr_cx = slide_pred_cancer)
df_calib$pr_ben = 1 - df_calib$pr_cx
## Add ISUP from clinical database (database_20190130.csv)
df_base <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'database_20190130.csv'))
# Create ISUP grade variable based on Gleason scores.    
df_base <- within(df_base, {
  ISUP <- NA_integer_
  ISUP[GS_all == '0 + 0'] <- 0
  ISUP[GS_all == '3 + 3'] <- 1
  ISUP[GS_all == '3 + 4'] <- 2
  ISUP[GS_all == '4 + 3'] <- 3
  ISUP[GS_all == '3 + 5'] <- 4
  ISUP[GS_all == '5 + 3'] <- 4
  ISUP[GS_all == '4 + 4'] <- 4
  ISUP[GS_all == '4 + 5'] <- 5
  ISUP[GS_all == '5 + 4'] <- 5
  ISUP[GS_all == '5 + 5'] <- 5
})
# remove duplicates
df_base <- df_base[!duplicated(df_base$slide), ]
# recode cx variable
# df_base$cx <- as.numeric(df_base$cx) - 1
df_base$cx <- ifelse(df_base$ISUP == 0, 0, 1)
df_base_calib <- filter(df_base, ext == '.svs')
df_base_test <- filter(df_base, ext == '.ndpi')
# Add variables to Test-Set
df_test <- left_join(df_test, select(df_base_test, c(slide, ISUP, cx, ext)), by = "slide")
# Add variables to Calib-Set
df_calib <- left_join(df_calib, select(df_base_calib, c(slide, ISUP, cx, ext)), by = "slide")

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
p_paired_aperio <- ks.test(df$X2, 'punif', 0, 1)$p.value

########################################################################################
## STHLM3-Paired: STHLM3-Aperio: Paired slides scanned on Hamamatsu
# Separate predictions on Calibrationset
########################################################################################
df_calib <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'df_test.csv'))
df_test <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'predictions_slide_rescanned_hamamatsu.csv'))

## Data-mgmt
# Rename columns Test-Set
df_test <- rename(df_test, X0=slide_pred_isup_0, X1=slide_pred_isup_1,
                  X2=slide_pred_isup_2, X3=slide_pred_isup_3,
                  X4=slide_pred_isup_4, X5=slide_pred_isup_5,
                  pr_cx = slide_pred_cancer)
df_test$pr_ben = 1 - df_test$pr_cx
# Rename columns Calib-Set
df_calib <- rename(df_calib, X0=slide_pred_isup_0, X1=slide_pred_isup_1,
                   X2=slide_pred_isup_2, X3=slide_pred_isup_3,
                   X4=slide_pred_isup_4, X5=slide_pred_isup_5,
                   pr_cx = slide_pred_cancer)
df_calib$pr_ben = 1 - df_calib$pr_cx
## Add ISUP from clinical database (database_20190130.csv)
df_base <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'database_20190130.csv'))
# Create ISUP grade variable based on Gleason scores.    
df_base <- within(df_base, {
  ISUP <- NA_integer_
  ISUP[GS_all == '0 + 0'] <- 0
  ISUP[GS_all == '3 + 3'] <- 1
  ISUP[GS_all == '3 + 4'] <- 2
  ISUP[GS_all == '4 + 3'] <- 3
  ISUP[GS_all == '3 + 5'] <- 4
  ISUP[GS_all == '5 + 3'] <- 4
  ISUP[GS_all == '4 + 4'] <- 4
  ISUP[GS_all == '4 + 5'] <- 5
  ISUP[GS_all == '5 + 4'] <- 5
  ISUP[GS_all == '5 + 5'] <- 5
})
# remove duplicates
df_base <- df_base[!duplicated(df_base$slide), ]
# recode cx variable
# df_base$cx <- as.numeric(df_base$cx) - 1
df_base$cx <- ifelse(df_base$ISUP == 0, 0, 1)
df_base_calib <- filter(df_base, ext == '.svs')
df_base_test <- filter(df_base, ext == '.ndpi')
# Add variables to Test-Set
df_test <- left_join(df_test, select(df_base_test, c(slide, ISUP, cx, ext)), by = "slide")
# Add variables to Calib-Set
df_calib <- left_join(df_calib, select(df_base_calib, c(slide, ISUP, cx, ext)), by = "slide")

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
p_paired_hamamatsu <- ks.test(df$X2, 'punif', 0, 1)$p.value

########################################################################################
## STHLM3-Paired: Evalauted on StGoran-slides (Aperio but different pathology lab)
# Separate predictions on Calibrationset
########################################################################################
df_calib <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'df_test.csv'))
df_test <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'predictions_slide_stgoran.csv'))

## Data-mgmt
# Rename columns Test-Set
df_test <- rename(df_test, X0=slide_pred_isup_0, X1=slide_pred_isup_1,
                  X2=slide_pred_isup_2, X3=slide_pred_isup_3,
                  X4=slide_pred_isup_4, X5=slide_pred_isup_5,
                  pr_cx = slide_pred_cancer)
df_test$pr_ben = 1 - df_test$pr_cx
# Rename columns Calib-Set
df_calib <- rename(df_calib, X0=slide_pred_isup_0, X1=slide_pred_isup_1,
                   X2=slide_pred_isup_2, X3=slide_pred_isup_3,
                   X4=slide_pred_isup_4, X5=slide_pred_isup_5,
                   pr_cx = slide_pred_cancer)
df_calib$pr_ben = 1 - df_calib$pr_cx
## Add ISUP from clinical database (database_20190130.csv)
df_base <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'database_20190130.csv'))
# Create ISUP grade variable based on Gleason scores.    
df_base <- within(df_base, {
  ISUP <- NA_integer_
  ISUP[GS_all == '0 + 0'] <- 0
  ISUP[GS_all == '3 + 3'] <- 1
  ISUP[GS_all == '3 + 4'] <- 2
  ISUP[GS_all == '4 + 3'] <- 3
  ISUP[GS_all == '3 + 5'] <- 4
  ISUP[GS_all == '5 + 3'] <- 4
  ISUP[GS_all == '4 + 4'] <- 4
  ISUP[GS_all == '4 + 5'] <- 5
  ISUP[GS_all == '5 + 4'] <- 5
  ISUP[GS_all == '5 + 5'] <- 5
})
# remove NDPI first
df_base <- filter(df_base, ext == '.svs')
# remove duplicates
df_base <- df_base[!duplicated(df_base$slide), ]
# recode cx variable
# df_base$cx <- as.numeric(df_base$cx) - 1
df_base$cx <- ifelse(df_base$ISUP == 0, 0, 1)
# Add variables to Calib-Set
df_calib <- left_join(df_calib, select(df_base, c(slide, ISUP, cx, ext)), by = "slide")
# Add variables to Test-Set
df_test <- left_join(df_test, select(df_base, c(slide, ISUP, cx, ext)), by = "slide")
df_test <- filter(df_test, ext == '.svs')

# Select only Aperio for comparison 
df_calib <- filter(df_calib, ext == '.svs')
table(df_calib$ext); table(df_test$ext)
# Select only ISUP4-5
df_calib <- filter(df_calib, ISUP %in% c(4,5))
df_test <- filter(df_test, ISUP %in% c(4,5))

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
p_paired_stgoran <- ks.test(df$X2, 'punif', 0, 1)$p.value

########################################################################################
## Combine
########################################################################################
p_values <- rbind(p_test, p_imagebase, p_ks, p_lars, p_paired_aperio, p_paired_hamamatsu, p_paired_stgoran)
p_values <- data.frame(p_values)
p_values

# Save
if (write_output){
  write.table(p_values, file = file.path('Output', 'Tables', 'kolmogorov','p_values_kolmogorov.csv'))
}
