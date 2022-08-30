
rm(list=ls())

library("knitr")
library("dplyr")
library("ggplot2")
library("gridExtra")
source("functions.R")

## Original Data
df_test <- read.csv(file = file.path('Data', 'Demo', 'sthlm3_aperio_hamamatsu.csv'))

########################################################################################
## Set Parameters
########################################################################################
write_output <- FALSE 

########################################################################################
## Add STHLM3_APERIO Calibration set
########################################################################################
df_calib <- read.csv(file = file.path('Data', 'Demo', 'sthlm3_aperio_calib.csv'))

########################################################################################
## Only use Aperio slides of Test-Set as CalibrationSet
########################################################################################
df_calib <- filter(df_calib, ext == '.svs')
table(df_calib$ext); table(df_test$ext)

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
cx <- calib_plot(pValues_cx, testSet = testSet) + ggtitle("Hamamatsu: CX")
# Cross-tabulate Outcome CX V. CP-CX Or Outcome ISUP V. CP-ISUP
df_pred01 <- df_predict_region(pValues_cx, sigfLevel = .001, outcome = 'CX')
df_pred05 <- df_predict_region(pValues_cx, sigfLevel = .005, outcome = 'CX')
df_pred1 <- df_predict_region(pValues_cx, sigfLevel = .01, outcome = 'CX')
tab01 <- tab_predict_region(df_pred01)
tab05 <- tab_predict_region(df_pred05)
tab1 <- tab_predict_region(df_pred1)
if (write_output){
  write.table(tab01, file = file.path('Output', 'Tables', 'hamamatsu_CX.csv'))
  write.table(tab05, file = file.path('Output', 'Tables', 'hamamatsu_CX.csv'), append = TRUE)
  write.table(tab1, file = file.path('Output', 'Tables', 'hamamatsu_CX.csv'), append = TRUE)
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
isup <- calib_plot(pValues_isup, testSet = testSet) + ggtitle("Hamamatsu: ISUP")
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
  write.table(tab05, file = file.path('Output', 'Tables', 'hamamatsu_ISUP.csv'))
  write.table(tab10, file = file.path('Output', 'Tables', 'hamamatsu_ISUP.csv'), append = TRUE)
  write.table(tab20, file = file.path('Output', 'Tables', 'hamamatsu_ISUP.csv'), append = TRUE)
  write.table(tab33, file = file.path('Output', 'Tables', 'hamamatsu_ISUP.csv'), append = TRUE)
  # Cancer-detection plot
  cx <- cx + ggtitle("")
  # Increase font-size
  cx <- cx + theme(text = element_text(size = 17))
  ggsave(filename = file.path('Output', 'Figures', 'cx_hamamatsu.png'), cx, 
         width = 6, 
         height = 6)
}
tab05;tab10;tab20;tab33

################################## end of program #################################
