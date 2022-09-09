rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
source("functions.R")

########################################################################################
## Set Parameters
########################################################################################
normalize <- FALSE
write_output <- FALSE 
use_testsplit <- FALSE

########################################################################################
## Test
########################################################################################
# II. Add test-split from I.
df_test <- read.csv(file = file.path('Data', 'Demo', 'df_test.csv'))
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
## Imagebase
df_test <- read.csv(file = file.path('Data',  'Demo', 'df_imagebase.csv'))
# Calibration-set
df_calib <- read.csv(file = file.path('Data', 'Demo', 'df_test.csv'))

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
df_test <- read.csv(file = file.path('Data', 'Demo', 'df_ks.csv'))

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
## Rare-subtypes
########################################################################################
# df1
rare_subtypes1 <- read.csv(file = file.path('Data',  'Demo', 
                                            'df_lars_special_cases_argmax.csv'))
# df2
rare_subtypes2 <- read.csv(file = file.path('Data',  'Demo',
                                            'rare_subtypes_ks.csv'))
# select matching variables
keepvars <- c('pr_cx', 'pr_ben', "cl_slide_cx", 'ISUP', 'cx', 'subtype')
rare_subtypes1 <- select(rare_subtypes1, keepvars)
rare_subtypes2 <- select(rare_subtypes2, keepvars)
df_test <- rbind(rare_subtypes1, rare_subtypes2)

# Calibration-set
df_calib <- read.csv(file = file.path('Data', 'Demo', 'df_test.csv'))

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
p_rare_subtypes <- ks.test(df$X2, 'punif', 0, 1)$p.value

########################################################################################
## STHLM3-Paired: STHLM3-Aperio: Paired slides scanned on Aperio
# Separate predictions on Calibrationset
########################################################################################
## Original Data
df_test <- read.csv(file = file.path('Data', 'Demo', 'sthlm3_aperio.csv'))
## Add STHLM3_APERIO Calibration set
df_calib <- read.csv(file = file.path('Data', 'Demo', 'sthlm3_aperio_calib.csv'))
df_calib <- filter(df_calib, ext == '.svs')

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
df_test <- read.csv(file = file.path('Data', 'Demo', 'sthlm3_aperio_hamamatsu.csv'))
## Add STHLM3_APERIO Calibration set
df_calib <- read.csv(file = file.path('Data', 'Demo', 'sthlm3_aperio_calib.csv'))
df_calib <- filter(df_calib, ext == '.svs')

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
## Stavanger - Test set 5
########################################################################################
df_test <- read.csv(file = file.path('Data', 'Demo', 'stavanger.csv'))
## Add STHLM3_APERIO Calibration set
df_calib <- read.csv(file = file.path('Data', 'Demo', 'sthlm3_aperio_calib.csv'))
# df_calib <- filter(df_calib, ext == '.svs')

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
p_stavanger <- ks.test(df$X2, 'punif', 0, 1)$p.value


########################################################################################
## Combine
########################################################################################
p_values <- rbind(p_test, p_imagebase, p_ks, p_rare_subtypes, p_paired_aperio, 
                  p_paired_hamamatsu, p_stavanger)
p_values <- data.frame(p_values)
p_values

# Save
if (write_output){
  write.table(p_values, file = file.path('Output', 'Tables', 'p_values_kolmogorov.csv'))
}
