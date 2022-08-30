#############################################################
### Common: This file contains the common functions which can
###         be shared among TCP and ICP. 

## Contains modified functions from the conformalClassification R-package
# https://cran.r-project.org/web/packages/conformalClassification/index.html
# Author: Niharika Gauraha and Ola Spjuth

## Additional functions for creation of prediction-regions / graphs / data management etc.
# Author Henrik Olsson
#############################################################

########################################################################################
# Computes efficiency of a conformal predictor, which is defined as the
# ratio of predictions with more than one class over the size of the testset
########################################################################################
CPEfficiency = function(matPValues, testLabels, sigfLevel = 0.05)
{
  if (is.null(matPValues) || is.null(testLabels)){
    stop("\n 'matPValues' and 'testLabels' are required as input\n")
  }
  
  nrTestCases = length(testLabels) #size of the test set
  
  signifTest = (matPValues  > sigfLevel)*1 #compute the prediction region
  
  err = 0
  for(i in 1:nrTestCases)
  {
    err = err + ( (sum(signifTest[i, ]) > 1) * 1 )
  }
  result = err/nrTestCases
  
  return(result)
}

########################################################################################
# Computes error rate of a conformal predictor, which is defined as
# the ratio of predictions with missing true class lables over the size of the testset
########################################################################################
CPErrorRate = function(matPValues, testLabels, sigfLevel = 0.05)
{
  if (is.null(matPValues) || is.null(testLabels)){
    stop("\n 'matPValues' and 'testLabels' are required as input\n")
  }
  
  nrTestCases = length(testLabels)
  
  signifTest = (matPValues  > sigfLevel)*1
  
  err = 0
  for(i in 1:nrTestCases)
  {
    err = err + (signifTest[i,testLabels[i]] == 0)*1
  }
  result = err/nrTestCases
  
  return(result)
}

########################################################################################
## Compute conformity scores and P-values: helper function
########################################################################################
computeConformityScores = function(calibrationSet = NULL)
{
  if(is.null(calibrationSet))
  {
    stop("Error: calibrationSet is NULL")
  }
  
  #The first colum should be the class labels
  calibLabels = as.numeric(calibrationSet[, 1])
  
  # predProb = predict(modelFit, calibrationSet[, -1], type="prob")
  predProb <- calibrationSet[-1]
  
  nrLabels = ncol(predProb) # number of class labels
  
  MCListConfScores = list() #Mondrian Class wise List of conformity scores
  for(i in 1:nrLabels)
  {
    classMembers = which(calibLabels == i)
    MCListConfScores[[i]] =  predProb[classMembers, i]
  }
  
  return(MCListConfScores)
}

computePValues = function(MCListConfScores = NULL, testConfScores = NULL)
{
  if(is.null(MCListConfScores) || is.null(testConfScores))
  {
    stop("Error: the list MCListConfScores is NULL")
    return(NULL)
  }
  
  nrTestCases = nrow(testConfScores)
  nrLabels = ncol(testConfScores)
  pValues = matrix(0, nrTestCases,  nrLabels)
  
  for(k in 1:nrTestCases)
  {
    for(l in 1:nrLabels)
    {
      alpha = testConfScores[k, l]
      classConfScores = MCListConfScores[[l]]
      pVal = length(which(classConfScores < alpha)) + runif(1)*length(which(classConfScores == alpha))
      pValues[k, l] = pVal/(length(classConfScores)+1)
    }
  }
  return(pValues)
}

ICPClassification = function(testSet, calibSet)
{
  MCListConfScores = computeConformityScores(calibSet)
  # testConfScores = predict(modelFit, testSet[, -1], type = "prob")
  testConfScores = testSet[-1]
  pValues = computePValues(MCListConfScores, testConfScores)
  
  return(pValues)
}

########################################################################################
## Evaluation: Conformal Prediction set [Single/Both/Null] at a given significance level 
########################################################################################
classmember <- function(pValues, sigfLevel = 0.05) 
{
  len <- dim(pValues)[1]
  signifTest <- pValues > sigfLevel
  n_empty <- sum(rowSums(signifTest) == 0)
  n_single <- sum(rowSums(signifTest) == 1)
  n_multiple <- sum(rowSums(signifTest) >= 2)
  
  return(list(N = c("empty" = n_empty, "single" = n_single, "multiple" = n_multiple), 
              Prop = c("empty" = n_empty/len, "single" = n_single/len, "multiple" = n_multiple/len)))
}

########################################################################################
## Stacked Area Graph: 
########################################################################################
area <- function(pval, lower=0, upper){
  # Empty/Single/Multiple by significance
  output <- matrix(0, 11, 3)
  for (i in 0:10) {
    prop <- classmember(pValues, sigfLevel = i*upper*.1)
    output[i+1,] <- prop$Prop
  }
  # Specificy plot parameters
  sig <- seq(lower, upper, upper*.1)
  prop <- as.vector(output)
  group <- factor(rep(seq(1:3), each = 11), label = c("Empty", "Single", "Multiple"))
  data <- data.frame(sig, prop, group)
  # Plot
  ggplot(data, aes(x=sig, y=prop, fill=group)) + 
    geom_area(alpha=0.6 , size=1, colour="black") + 
    scale_y_continuous(breaks = c(seq(0, 1, .1))) + 
    scale_x_continuous(breaks = c(seq(lower, upper, upper*.1))) +
    xlab("significance") + 
    ylab("Proportion")
}

########################################################################################
## Calibration plot (using ggplot, easier to combine plots than base R)
########################################################################################
calib_plot <- function(pValues, testSet, color="blue")
{
  if (is.null(pValues) || is.null(testSet)){
    stop("\n 'pValues' and 'testSet' are required as input\n")
  }
  
  testLabels = testSet[, 1] #True class labels
  nrTestCases = length(testLabels)
  
  #compute error rate for the range of significance levels
  sigLevels = seq(0,1, .01)
  errorRate = rep(1, length(sigLevels))
  for(i in 1:(length(sigLevels) - 1)){
    signifTest = (pValues  > sigLevels[i])*1
    err = 0
    for(j in 1:nrTestCases)
    {
      err = err + ( (signifTest[j, testLabels[j]] == 0)*1 )
    }
    
    errorRate[i] = err/nrTestCases
  }
  
  data <- data.frame(sigLevels, errorRate)
  p <- ggplot(data, aes(x=sigLevels, y=errorRate)) + 
    geom_line() + 
    geom_abline(colour = "red") + 
    scale_x_continuous(breaks = seq(0, 1, 0.2)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    xlab("Significance") +
    ylab("Error rate")
  
  return(p)
}

########################################################################################
## Data-Mgmt
########################################################################################
generate_vars <- function(outcome, testdata){
  if (outcome == "CX") {
    ## make sure first column is always the label and class labels are always 1, 2, ...
    # Calibration Set
    df_calib <- select(df_calib, c(cx, pr_ben, pr_cx))
    df_calib$cx <- ifelse(df_calib$cx == 0, 1, 2)
    # Test Set
    df_test <- select(df_test, c(cx, pr_ben, pr_cx))
    df_test$cx <- ifelse(df_test$cx == 0, 1, 2)
    col_outcome <- c(1, 2)
  } 
  if (outcome == "ISUP") {
    ## make sure first column is always the label and class labels are always 1, 2, ...
    # Calibration Set
    df_calib <- select(df_calib, c(ISUP, X0:X5))
    df_calib$ISUP <- df_calib$ISUP+1
    # Test Set
    df_test <- select(df_test, c(ISUP, X0:X5))
    df_test$ISUP <- df_test$ISUP+1
    col_outcome <- c(1, 2, 3, 4, 5, 6)
  }
  if (outcome == "ISUP_group") {
    ## make sure first column is always the label and class labels are always 1, 2, ...
    # Calibration Set
    df_calib <- select(df_calib, c(ISUP, X0:X3))
    df_calib$ISUP <- df_calib$ISUP+1
    # Test Set
    df_test <- select(df_test, c(ISUP, X0:X3))
    df_test$ISUP <- df_test$ISUP+1
    col_outcome <- c(1, 2, 3, 4)
  }
  
  if (testdata == "test"){
    testSet = df_test
    testLabels <- df_test[, 1]
  }
  if (testdata == "ks"){
    testSet = df_ks1
    testLabels <- df_ks1[, 1]
  }
  if (testdata == "ks_cgan"){
    testSet = df_KS1_CGAN
    testLabels <- df_KS1_CGAN[, 1]
    df_calib <- df_calib
  }
  return(list(testSet = testSet, testLabels = testLabels, df_calib = df_calib, col_outcome = col_outcome))
}

########################################################################################
## Create data with Prediction-Region
########################################################################################
df_predict_region <- function(pValues, sigfLevel, outcome){
  # Predicted classes
  signifTest <- pValues > sigfLevel
  # Number of predicted labels in each prediction set
  sum <- rowSums(signifTest)
  # Error predictions: True label missing in prediction set
  error <- c()
  for (i in 1:length(testLabels)){
    error[i] <- signifTest[i, testLabels[i]] == 0
  }
  # Prediction regions: Significant Labels
  signifTest <- as.data.frame(signifTest)
  names(signifTest) <- col_outcome
  signifTest$true_cols <- apply(signifTest, 1, function(data)
    names(which(data == T)))
  signifTest$true_cols <- as.character(signifTest$true_cols)
  signifTest$true_cols[signifTest$true_cols == "character(0)"] <- "Empty"
  
  # combine
  df <- data.frame(testLabels, signifTest, sum, error)
  
  # Create Empty / Error / Single / Multiple
  df <- within(df, {
    pred_group <- NA
    pred_group[sum == 0] <- 1
    pred_group[sum != 0 & error == TRUE] <- 2
    pred_group[sum == 1 & error == FALSE] <- 3
    pred_group[sum > 1 & error == FALSE] <- 4
    error <- ifelse(sum == 0, FALSE, error)
    
    empty <- NA
    empty <- ifelse(sum == 0, TRUE, FALSE)
    
    if (outcome == 'ISUP'){
      testLabels <- factor(testLabels,  levels = c(1,2,3,4,5,6), 
                           labels = c("Benign", "ISUP 1", "ISUP 2", "ISUP 3", "ISUP 4", "ISUP 5"))
      pred_group <- factor(pred_group, levels = c(1,2,3,4), labels = c("Empty", "Error", "Single", "Multiple"))
    }
    if (outcome == 'ISUP_group'){
      testLabels <- factor(testLabels,  levels = c(1,2,3,4), 
                           labels = c("Benign", "ISUP 1", "ISUP 2-3", "ISUP 4-5"))
      pred_group <- factor(pred_group, levels = c(1,2,3,4), labels = c("Empty", "Error", "Single", "Multiple"))
    }
    if (outcome == 'CX'){
      testLabels <- factor(testLabels, levels = c(1,2), labels = c("Benign", "CX"))
      pred_group <- factor(pred_group, levels = c(1,2,3, 4), labels = c("Empty", "Error", "Single", "Multiple"))
    }
  })
  return(df)
}

########################################################################################
## Output Table: Prediction-Set by ISUP-grade
########################################################################################
tab_predict_region <- function(df){
  x <- df
  tab <- table(x$pred_group, x$testLabels, useNA="ifany")
  tab_margins <- addmargins(tab, 2)
  mytable <- addmargins(tab)
  for (i in 1: dim(mytable)[1]-1){
    mytable[i,] <- paste(mytable[i,], " (", round( prop.table(tab_margins,2)[i,],2), ")", sep="")
  }
  return(mytable)
}

########################################################################################
## Output: Table of prediction regions by clinical diagnosis
# Input: df with prediction regions and clinical diagnosis
########################################################################################
library("stringr")

predgroups <- function(df_pred){
  x <- df_pred
  # Filter out multiple predictions and recode ISUP
  x <- filter(x, pred_group == "Multiple")
  x$true_cols <- gsub("[^0-9 ]", "", x$true_cols)
  x$true_cols <- str_replace(x$true_cols, "2", "ISUP1")
  x$true_cols <- str_replace(x$true_cols, "3", "ISUP2")
  x$true_cols <- str_replace(x$true_cols, "4", "ISUP3")
  x$true_cols <- str_replace(x$true_cols, "5", "ISUP4")
  x$true_cols <- str_replace(x$true_cols, "6", "ISUP5")

  # Tabulate prediction-regions by ISUP
  tab <- table(x$true_cols, x$testLabels, useNA="ifany")
  tab_margins <- addmargins(tab, 2)
  mytable <- addmargins(tab)
  for (i in 1: dim(mytable)[1]-1){
    mytable[i,] <- paste(mytable[i,], " (", round( prop.table(tab_margins,2)[i,],2), ")", sep="")
  }
  return(mytable)
}

