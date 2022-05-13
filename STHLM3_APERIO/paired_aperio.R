
rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("pROC")
source("Programs/functions.R")

## Original Data
df_ks1 <- read.csv(file = file.path('Data', 'Derived', 'conformal_20200225', 'df_KS1.csv'))
df_KS1_CGAN <- read.csv(file = file.path('Data', 'Derived', 'conformal_20200225', 'df_KS1_CGAN.csv'))
# df_test <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'models_20201202', 'predictions_slide_test.csv'))
df_test <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'predictions_slide_rescanned_aperio.csv'))

########################################################################################
## Set Parameters
########################################################################################
normalize <- FALSE
ratioTrain <- 0.5
create_testsplit <- FALSE
write_output <- FALSE 

########################################################################################
## Sample part of Internal Test-Data as calibration
########################################################################################
# I. Run only once
if (create_testsplit){
  df_calib <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'df_test.csv'))
  result <- sample(1:nrow(df_calib), ratioTrain*nrow(df_calib))
  df_calib$calibset <- ifelse(as.numeric(rownames(df_calib)) %in% result, 0, 1)
  write.csv(df_calib, file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'df_test.csv'))
}

# II. Add test-split from I.
df_calib <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'df_test.csv'))
# df_calib <- filter(df_calib, calibset == 1)

########################################################################################
## Data-mgmt
# Rename columns to match code and functions
# Add ISUP and CX from clinical database (database_20190130.csv)
########################################################################################
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

########################################################################################
## Only use Aperio slides of Test-Set as CalibrationSet

# If they are similar, might be an issue with poor training and loss in performance in 
  # DNN-training
########################################################################################
df_calib <- filter(df_calib, ext == '.svs')
table(df_calib$ext); table(df_test$ext)

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
# Create Empty / Error / Single / Multiple
df_calib$slide_class_isup <- factor(df_calib$slide_class_isup, levels = c(0, 1, 2, 3, 4, 5)); 
df_test$slide_class_isup <- factor(df_test$slide_class_isup, levels = c(0, 1, 2, 3, 4, 5))

caret::confusionMatrix(factor(df_calib$ISUP), df_calib$slide_class_isup); vcd::Kappa(table(df_calib$ISUP, df_calib$slide_class_isup))
caret::confusionMatrix(factor(df_test$ISUP), df_test$slide_class_isup); vcd::Kappa(table(df_test$ISUP, df_test$slide_class_isup))

with(df_calib, roc(cx, pr_cx)); with(df_test, roc(cx, pr_cx))

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
cx <- calib_plot(pValues_cx, testSet = testSet) + ggtitle("Aperio: CX")
# Cross-tabulate Outcome CX V. CP-CX Or Outcome ISUP V. CP-ISUP
df_pred01 <- df_predict_region(pValues_cx, sigfLevel = .001, outcome = 'CX')
df_pred05 <- df_predict_region(pValues_cx, sigfLevel = .005, outcome = 'CX')
df_pred1 <- df_predict_region(pValues_cx, sigfLevel = .01, outcome = 'CX')
tab01 <- tab_predict_region(df_pred01)
tab05 <- tab_predict_region(df_pred05)
tab1 <- tab_predict_region(df_pred1)
if (write_output){
  write.table(tab01, file = file.path('Output', 'Tables', 'PredictionSets', 'STHLM3_APERIO', 'aperio_CX.csv'))
  write.table(tab05, file = file.path('Output', 'Tables', 'PredictionSets', 'STHLM3_APERIO', 'aperio_CX.csv'), append = TRUE)
  write.table(tab1, file = file.path('Output', 'Tables', 'PredictionSets', 'STHLM3_APERIO', 'aperio_CX.csv'), append = TRUE)
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
isup <- calib_plot(pValues_isup, testSet = testSet) + ggtitle("Aperio: ISUP")
# Cross-tabulate Outcome CX V. CP-CX Or Outcome ISUP V. CP-ISUP
df_pred05 <- df_predict_region(pValues_isup, sigfLevel = .05, outcome = 'ISUP')
df_pred10 <- df_predict_region(pValues_isup, sigfLevel = .1, outcome = 'ISUP')
df_pred20 <- df_predict_region(pValues_isup, sigfLevel = .2, outcome = 'ISUP')
df_pred33 <- df_predict_region(pValues_isup, sigfLevel = 1/3, outcome = 'ISUP')
tab05 <- tab_predict_region(df_pred05)
tab10 <- tab_predict_region(df_pred10)
tab20 <- tab_predict_region(df_pred20)
tab33 <- tab_predict_region(df_pred33)
if (write_output){
  write.table(tab05, file = file.path('Output', 'Tables', 'PredictionSets', 'STHLM3_APERIO', 'aperio_ISUP.csv'))
  write.table(tab10, file = file.path('Output', 'Tables', 'PredictionSets', 'STHLM3_APERIO', 'aperio_ISUP.csv'), append = TRUE)
  write.table(tab20, file = file.path('Output', 'Tables', 'PredictionSets', 'STHLM3_APERIO', 'aperio_ISUP.csv'), append = TRUE)
  write.table(tab33, file = file.path('Output', 'Tables', 'PredictionSets', 'STHLM3_APERIO', 'aperio_ISUP.csv'), append = TRUE)
  # Cancer-detection plot
  # cx <- cx + ggtitle("Same scanner")
  # ggsave(filename = file.path('Output', 'Plots', 'calibration', 'cx_aperio.png'), cx, 
  #        width = 2, 
  #        height = 2)
  cx <- cx + ggtitle("")
  # Increase font-size
  cx <- cx + theme(text = element_text(size = 17))  
  ggsave(filename = file.path('Output', 'Plots', 'calibration', 'nolabels', 'cx_aperio.png'), cx, 
         width = 6, 
         height = 6)
}
tab05;tab10;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

########################################################################################
## Save Figure
########################################################################################
# calibration_aperio <- grid.arrange(cx, isup, ncol = 2)
# ggsave(calibration_aperio, file = file.path('Output', 'Plots', 'calibration_aperio.png'))

################################## end of program #################################
