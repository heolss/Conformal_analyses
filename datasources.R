rm(list=ls())

library("dplyr")

########################################################################################
## Set Parameters
########################################################################################
write_output <- FALSE 

########################################################################################
## Data sources
########################################################################################
# Train
df_train <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_train.csv'))
# CalibrationSet
df_test <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_test.csv'))
df_calib <- filter(df_test, calibset == 1)
# TestSet
df_test <- filter(df_test, calibset == 0)
## Imagebase
df_imagebase <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_imagebase.csv'))
# Imagebase: Add MODE as ISUP
imagebase <- read.csv(file = file.path('Data', 'Derived', 'CP_LO', 'FinalResultsProstate_with_AI_anonymized.csv'))
df_imagebase <- left_join(df_imagebase, select(imagebase, c(Path.no, MODE, Consensus)), by = c("slide" = "Path.no"))
df_imagebase <- filter(df_imagebase, !is.na(MODE))
df_imagebase$ISUP <- df_imagebase$MODE
## Lars-special-cases
df_LSC <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_lars_special_cases_argmax.csv'))
df_LSC <- read.csv(file = file.path('Data', 'Derived', 'CP_LO_argmaxvoting', 'df_lars_special_cases_argmax.csv'))
# Recode error in clinical data
# OBS! These needs to be checked with Lars
table(df_LSC$cx); table(df_LSC$ISUP)
df_LSC$ISUP[df_LSC$ISUP == 5] <- 1
df_LSC$ISUP[df_LSC$slide == 'SP845_2012 1'] <- 1
df_LSC$ISUP[df_LSC$slide == 'SP1475_2011 6'] <- 1
df_LSC$ISUP[df_LSC$slide == 'SP1034_2012 1V'] <- 1
# df_LSC$ISUP[df_LSC$slide == 'SP960_2013 1A'] <- 1
df_LSC$comment <- as.character(df_LSC$comment)
df_LSC$comment[df_LSC$slide == 'SP960_2013 1A'] <- "PAH"
df_LSC$cx <- ifelse(df_LSC$ISUP >= 1, 1, 0)
table(df_LSC$cx); table(df_LSC$ISUP)
# Remove RP 
df_LSC <- filter(df_LSC, slide != 'SP1034_2012 1V')
# KS
df_KS <- read.csv(file = file.path('Data', 'Derived', 'CP_LO', 'df_ks.csv'))
########################################################################################
## STG
## Data-mgmt
# Rename columns to match code and functions
# Add ISUP and CX from clinical database (database_20190130.csv)
########################################################################################
df_STG <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'predictions_slide_stgoran.csv'))
# Rename columns Test-Set
df_STG <- rename(df_STG, X0=slide_pred_isup_0, X1=slide_pred_isup_1,
                 X2=slide_pred_isup_2, X3=slide_pred_isup_3,
                 X4=slide_pred_isup_4, X5=slide_pred_isup_5,
                 pr_cx = slide_pred_cancer)
df_STG$pr_ben = 1 - df_STG$pr_cx

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
df_base$cx <- ifelse(df_base$ISUP == 0, 0, 1)

# Add variables to Test-Set
df_STG <- left_join(df_STG, select(df_base, c(slide, subject, ISUP, cx, ext)), by = "slide")
df_STG <- filter(df_STG, ext == '.svs')
df_STG <- filter(df_STG, ISUP %in% c(4,5))
df_STG <- rename(df_STG, man = subject)

########################################################################################
## Paired-Aperio-Hamamatsu
## Data-mgmt
# Rename columns to match code and functions
# Add ISUP and CX from clinical database (database_20190130.csv)
########################################################################################
df_paired <- read.csv(file = file.path('Data', 'Derived', 'STHLM3_APERIO', 'predictions_slide_rescanned_hamamatsu.csv'))
# Data-mgmt
# Rename columns Test-Set
df_paired <- rename(df_paired, X0=slide_pred_isup_0, X1=slide_pred_isup_1,
                    X2=slide_pred_isup_2, X3=slide_pred_isup_3,
                    X4=slide_pred_isup_4, X5=slide_pred_isup_5,
                    pr_cx = slide_pred_cancer)
df_paired$pr_ben = 1 - df_paired$pr_cx
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
df_base_test <- filter(df_base, ext == '.ndpi')
# Add variables to Test-Set
df_paired <- left_join(df_paired, select(df_base_test, c(slide, subject, ISUP, cx, ext)), by = "slide")
df_paired <- rename(df_paired, man = subject)

########################################################################################
## Combine
########################################################################################
df_train$source <- 1
df_calib$source <- 2
df_test$source <- 3
df_imagebase$source <- 4
df_KS$source <- 5
df_LSC$source <- 6
df_STG$source <- 7
df_paired$source <- 8

df <- rbind(select(df_train, c(slide, man, cx, ISUP, source)),
            select(df_calib, c(slide, man, cx, ISUP, source)),
            select(df_test, c(slide, man, cx, ISUP, source)),
            select(df_imagebase, c(slide, man, cx, ISUP, source)),
            select(df_KS, c(slide, man, cx, ISUP, source)),
            select(df_LSC, c(slide, man, cx, ISUP, source)),
            select(df_STG, c(slide, man, cx, ISUP, source)),
            select(df_paired, c(slide, man, cx, ISUP, source)))

df$source <- factor(df$source, levels = c(1,2,3,4,5,6,7,8),
                    labels = c("Train", "CalibrationSet", "TestSet", "Imagebase", "KS", "Rare prostate tissue morphology", "STG", "Paired"))

df$ISUP <- factor(df$ISUP,  levels = c(0,1,2,3,4,5), 
                  labels = c("Benign", "ISUP 1", "ISUP 2", "ISUP 3", "ISUP 4", "ISUP 5"))

########################################################################################
## Table: Tabulate data source by ISUP
########################################################################################
tab <- function(df){
  x <- df
  tab <- table(x$ISUP, x$source, useNA="ifany")
  tab_margins <- addmargins(tab, 2)
  mytable <- addmargins(tab)
  for (i in 1: dim(mytable)[1]-1){
    mytable[i,] <- paste(mytable[i,], " (", round( prop.table(tab_margins,2)[i,],2), ")", sep="")
  }
  return(mytable)
}

datasources <- tab(df)
datasources


if (write_output){
  write.table(datasources, file = file.path('Output', 'Tables', 'Datasources', 'datasources.csv'))
}

################################## end of script #################################

########################################################################################
## Check overlap
########################################################################################
# a <- filter(df_LSC, slide %in% df_train$slide)
# b <- filter(df_train, slide %in% df_LSC$slide)
# 
# c <- filter(df_imagebase, slide %in% df_calib$slide)
# 
# dim(filter(df_imagebase, slide %in% df_calib$slide))[1]
# dim(filter(df_imagebase, slide %in% df_test$slide))[1]
