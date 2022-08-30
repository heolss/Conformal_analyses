rm(list=ls())

library("dplyr")

########################################################################################
## Set Parameters
########################################################################################
write_output <- FALSE 

########################################################################################
## Data sources
########################################################################################
## Train
df_train <- read.csv(file = file.path('Data',  'Demo', 'df_train.csv'))
## CalibrationSet
df_test <- read.csv(file = file.path('Data',  'Demo', 'df_test.csv'))
df_calib <- filter(df_test, calibset == 1)
## TestSet
df_test <- filter(df_test, calibset == 0)
## Imagebase
df_imagebase <- read.csv(file = file.path('Data',  'Demo', 'df_imagebase.csv'))
## Rare-subtypes
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
df_rare_subtypes <- rbind(rare_subtypes1, rare_subtypes2)
## KS
df_KS <- read.csv(file = file.path('Data',  'Demo', 'df_ks.csv'))
## Paired-Aperio-Hamamatsu
df_paired <- read.csv(file = file.path('Data',  'Demo', 'sthlm3_aperio_hamamatsu.csv'))

########################################################################################
## Combine
########################################################################################
df_train$source <- 1
df_calib$source <- 2
df_test$source <- 3
df_imagebase$source <- 4
df_KS$source <- 5
df_rare_subtypes$source <- 6
df_paired$source <- 7

df <- rbind(select(df_train, c(cx, ISUP, source)),
            select(df_calib, c(cx, ISUP, source)),
            select(df_test, c(cx, ISUP, source)),
            select(df_imagebase, c(cx, ISUP, source)),
            select(df_KS, c(cx, ISUP, source)),
            select(df_rare_subtypes, c(cx, ISUP, source)),
            select(df_paired, c(cx, ISUP, source)))

df$source <- factor(df$source, levels = c(1,2,3,4,5,6,7),
                    labels = c("Train", "CalibrationSet", "TestSet", "Imagebase", "KS", 
                               "Rare prostate tissue morphology", "Paired"))

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
  write.table(datasources, file = file.path('Output', 'Tables', 'datasources.csv'))
}

################################## end of script #################################

