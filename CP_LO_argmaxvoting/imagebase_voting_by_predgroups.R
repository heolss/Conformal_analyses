########################################################################################
# Make table of prediction regions for different confidence levels - additionally compute the coverage of 
# IMB-panel votes that are covered by prediction regions.

# Table No. votes / freq lying in our prediction-sets (below different alternatives)
# Prediction-regions V. IMB-votes: Report overlapping, size of sets

# Cumulative frequency of votes for the given prediction-set: 
# e.g. ISUP1 -> % of ISUP1 by IMB. ISUP1-ISUP2 -> cumulative-% of ISUP1-ISUP2
########################################################################################
rm(list=ls())

set.seed(12)

library("dplyr")
library("tidyr")
library("ggplot2")
library("gridExtra")
source("functions.R")
library("reshape2")

########################################################################################
## Run Conformal prediction on Imagebase data and recode variables

# Generate Table of predgroups by ISUP
# Input to # Combine with Imbase-panel below
########################################################################################
source("CP_LO_argmaxvoting/CP_LO_IMAGEBASE.R")
# DF_PRED <- df_pred20
DF_PRED <- df_pred33

########################################################################################
## Imagebase data
########################################################################################
imagebase <- read.csv(file = file.path('Data', 'Demo', 'df_imagebase.csv'))

########################################################################################
## Relative frequency of ISUP-group-votes by biopsy-core
########################################################################################
# rename ID-variable
imagebase <- plyr::rename(imagebase, c("Path.no" = "slide"))
# Each Imagebase-panel-member
vars <- c(grep("IB", names(imagebase), value = TRUE, ignore.case = TRUE))

## Relative frequency of ISUP-assessment of Imagebase-panel for each biopsy
# Melt by Patient
long <- melt(imagebase,
             id.vars = c("slide", "MODE"),
             measure.vars=c(vars),
             variable.name="name",
             value.name="value")
# Relative frequency by biopsy
long <- long %>%
  group_by(slide, MODE, value) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
long <- select(long, c(slide, value, freq))
# long to wide
wide <- dcast(long, slide + MODE ~ value, value.var="freq")
wide
df_imb <- wide

########################################################################################
# Combine with Imbase-panel
########################################################################################
# Add slide variable
df <- cbind(df_test$slide, DF_PRED)
df <- plyr::rename(df, c("df_test$slide" = "slide"))
# Select variables
df <- select(df, c(slide, true_cols, sum))
# Add IMB-predictions
df <- left_join(df, df_imb, by = "slide")

df$true_cols <- gsub("[^0-9 ]", "", df$true_cols)

########################################################################################
# Cleaon up ISUP-coding
########################################################################################
label <- function(data){
  x <- data
  x$true_cols <- gsub("[^0-9 ]", "", x$true_cols)
  x$true_cols <- str_replace(x$true_cols, "2", "ISUP1")
  x$true_cols <- str_replace(x$true_cols, "3", "ISUP2")
  x$true_cols <- str_replace(x$true_cols, "4", "ISUP3")
  x$true_cols <- str_replace(x$true_cols, "5", "ISUP4")
  x$true_cols <- str_replace(x$true_cols, "6", "ISUP5")
  
  x <- plyr::rename(x, c("1" = "ISUP1"))
  x <- plyr::rename(x, c("2" = "ISUP2"))
  x <- plyr::rename(x, c("3" = "ISUP3"))
  x <- plyr::rename(x, c("4" = "ISUP4"))
  x <- plyr::rename(x, c("5" = "ISUP5"))
  return(x)
}
df <- label(df)
head(df, 10)

########################################################################################
# Select rows that overlap: Our prediction-regions & IMB-votes: 
# Among overlapping rows: Calculate IMB-voting
# Calculate IMB-votes / frequency of votes for our prediction-sets
########################################################################################
long <- melt(df,
             id.vars = c("slide", "true_cols"),
             measure.vars=c("ISUP1","ISUP2","ISUP3","ISUP4","ISUP5"),
             variable.name="name",
             value.name="value")
# select only rows with values
long <- long %>% 
  filter(!is.na(value)) %>% 
  arrange(slide)
head(long, 10)
# Select overlapping rows
long$true_cols <- gsub(" ", "|", long$true_cols)
a <- filter(long, str_detect(name, true_cols))
# Summarize over slides
b <- a %>%
  group_by(slide) %>%
  summarise(Frequency = sum(value))
# combine
c <- left_join(df, b, by = "slide")
c$Frequency <- ifelse(is.na(c$Frequency), 0, c$Frequency)
summary(c$Frequency); IQR(c$Frequency)
hist(c$Frequency)

df <- c

########################################################################################
# 1. Table prediction-sets over MODE at 67% Confidence-level
# 2. Median IQR of all IMB-votes for prediction-set
########################################################################################
# Our predicionts
tab <- table(df$true_cols, df$MODE, useNA="ifany")
tab_margins <- addmargins(tab, 2)
mytable <- addmargins(tab)
for (i in 1: dim(mytable)[1]-1){
  mytable[i,] <- paste(mytable[i,], " (", round( prop.table(tab_margins,2)[i,],2), ")", sep="")
}
t1 <- mytable
t2 <- df %>%
  group_by(true_cols) %>%
  summarise(iqr = round(IQR(Frequency), 2), Frequency = round(median(Frequency), 2))
t1; t2

# Overall coverage - frequency of votes for our prediction-sets
coverage <- summary(df$Frequency)['Median'][[1]]

if (write_output){
  write.table(t1, file = file.path('Output', 'Tables', 'predgroup_imb_voting.csv'))
  write.table(t2, file = file.path('Output', 'Tables', 'predgroup_imb_voting.csv'), append = TRUE)
  write.table(coverage, file = file.path('Output', 'Tables', 'predgroup_imb_voting.csv'), append = TRUE)
}

################################## end of program #################################
