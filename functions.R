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
# Computes observed fuzziness, which is defined as
# the sum of all p-values for the incorrect class labels.
########################################################################################
CPObsFuzziness = function(matPValues, testLabels)
{
  if (is.null(matPValues) || is.null(testLabels)){
    stop("\n 'matPValues' and 'testLabels' are required as input\n")
  }
  
  nrTestCases = length(testLabels)
  
  sumPValues = 0
  for(indxTestSet in 1:nrTestCases)
  {
    exclude = testLabels[indxTestSet] #exclude the p-value of the true label
    sumPValues = sumPValues + sum(matPValues[indxTestSet, -exclude])
  }
  result = sumPValues/nrTestCases
  return(result)
}

########################################################################################
# Computes the deviation from exact validity as the Euclidean norm of
# the difference of the observed error and the expected error
########################################################################################
CPValidity = function(matPValues = NULL, testLabels = NULL)
{
  if (is.null(matPValues) || is.null(testLabels)){
    stop("\n 'matPValues' and 'testLabels' are required as input\n")
  }
  
  signifSet = seq(.01, .99, by=.01) #significance level set
  
  nrTestCases = length(testLabels)
  errAtSignif = rep(0, length(signifSet))
  
  for(indx in 1: length(signifSet)){
    signifTest = (matPValues  > signifSet[indx])*1
    
    err = 0
    for(i in 1:nrTestCases)
    {
      err = err + ( (signifTest[i, testLabels[i]] == 0) * 1 )
    }
    err = err/nrTestCases
    errAtSignif[indx] = (err - signifSet[indx])^2
  }
  
  result = sqrt(sum(errAtSignif))
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
  
  MCListConfScores = list() #Moderian Class wise List of conformity scores
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
    # Test Set
    df_ks1 <- select(df_ks1, c(cx, pr_ben, pr_cx))
    df_ks1$cx <- ifelse(df_ks1$cx == 0, 1, 2)
    # CGAN
    df_KS1_CGAN <- select(df_KS1_CGAN, c(cx, pr_ben, pr_cx))
    df_KS1_CGAN$cx <- ifelse(df_KS1_CGAN$cx == 0, 1, 2)
    # Outcome-Test-classes
    col_outcome <- c(1, 2)
    # Calibrated
    # df_calib_calibrated <- select(df_calib_calibrated, c(cx, pr_ben, pr_cx))
    # df_calib_calibrated$cx <- ifelse(df_calib_calibrated$cx == 0, 1, 2)
    # df_ks1_calibrated <- select(df_ks1_calibrated, c(cx, pr_ben, pr_cx))
    # df_ks1_calibrated$cx <- ifelse(df_ks1_calibrated$cx == 0, 1, 2)
    # df_KS1_CGAN_calibrated <- select(df_KS1_CGAN_calibrated, c(cx, pr_ben, pr_cx))
    # df_KS1_CGAN_calibrated$cx <- ifelse(df_KS1_CGAN_calibrated$cx == 0, 1, 2)
  } 
  if (outcome == "ISUP") {
    ## make sure first column is always the label and class labels are always 1, 2, ...
    # Calibration Set
    df_calib <- select(df_calib, c(ISUP, X0:X5))
    df_calib$ISUP <- df_calib$ISUP+1
    # Test Set
    df_test <- select(df_test, c(ISUP, X0:X5))
    df_test$ISUP <- df_test$ISUP+1
    # KS-Data
    df_ks1 <- select(df_ks1, c(ISUP, X0:X5))
    df_ks1$ISUP <- df_ks1$ISUP+1
    # CGAN
    df_KS1_CGAN <- select(df_KS1_CGAN, c(ISUP, X0:X5))
    df_KS1_CGAN$ISUP <- df_KS1_CGAN$ISUP+1
    # Outcome-Test-classes
    col_outcome <- c(1, 2, 3, 4, 5, 6)
    # Calibrated
    # df_calib_calibrated <- select(df_calib_calibrated, c(ISUP, X0:X5))
    # df_calib_calibrated$ISUP <- df_calib_calibrated$ISUP+1
    # df_ks1_calibrated <- select(df_ks1_calibrated, c(ISUP, X0:X5))
    # df_ks1_calibrated$ISUP <- df_ks1_calibrated$ISUP+1
    # df_KS1_CGAN_calibrated <- select(df_KS1_CGAN_calibrated, c(ISUP, X0:X5))
    # df_KS1_CGAN_calibrated$ISUP <- df_KS1_CGAN_calibrated$ISUP+1
  }
  if (outcome == "ISUP_group") {
    ## make sure first column is always the label and class labels are always 1, 2, ...
    # Calibration Set
    df_calib <- select(df_calib, c(ISUP, X0:X3))
    df_calib$ISUP <- df_calib$ISUP+1
    # Test Set
    df_test <- select(df_test, c(ISUP, X0:X3))
    df_test$ISUP <- df_test$ISUP+1
    # KS-Data
    df_ks1 <- select(df_ks1, c(ISUP, X0:X3))
    df_ks1$ISUP <- df_ks1$ISUP+1
    # CGAN
    df_KS1_CGAN <- select(df_KS1_CGAN, c(ISUP, X0:X3))
    df_KS1_CGAN$ISUP <- df_KS1_CGAN$ISUP+1
    # Outcome-Test-classes
    col_outcome <- c(1, 2, 3, 4)
    # Calibrated
  }
  
  if (testdata == "test"){
    testSet = df_test
    testLabels <- df_test[, 1]
  }
  if (testdata == "ks"){
    testSet = df_ks1
    testLabels <- df_ks1[, 1]
  }
  # if (testdata == "ks_calibrated"){
  #   testSet = df_ks1_calibrated
  #   testLabels <- df_ks1_calibrated[, 1]
  #   df_calib <- df_calib_calibrated
  # }
  if (testdata == "ks_cgan"){
    testSet = df_KS1_CGAN
    testLabels <- df_KS1_CGAN[, 1]
    df_calib <- df_calib
  }
  # if (testdata == "ks_cgan_calibrated"){
  #   testSet = df_KS1_CGAN_calibrated
  #   testLabels <- df_KS1_CGAN_calibrated[, 1]
  #   df_calib <- df_calib_calibrated
  # }
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
    # pred_group <- factor(pred_group, labels = c("Empty", "Error", "Single", "Multiple"))
    # Recode Empty predictions from error predictions
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
## Logistic calibration Function
########################################################################################
# OBS! It is possible that calibration should be performed only on k-1 groups? Remove p0?
logistic_calibration <- function(data, outcome){
  if (outcome == "CX"){
    # lp1 <- with(data, log(pr_cx / (1 - pr_cx)))
    # fit1 <- glm(cx == 1 ~ lp1, family = binomial, data = data)
    # p1 <- predict(fit1, type = "response")
    # # Updated Predictions with calibrated probabilities
    # data <- data.frame(data, p1)
    # # Normalize predictions:
    # data <- data %>%
    #   mutate(pr_cx = p1,
    #          pr_ben = 1 - p1) %>%
    #   select(-c(p1))
    lp0 <- with(data, log(pr_ben / (1 - pr_ben)))
    fit0 <- glm(cx == 0 ~ lp0, family = binomial, data = data)
    p0 <- predict(fit0, type = "response")
    lp1 <- with(data, log(pr_cx / (1 - pr_cx)))
    fit1 <- glm(cx == 1 ~ lp1, family = binomial, data = data)
    p1 <- predict(fit1, type = "response")

    # Updated Predictions with calibrated probabilities
    data <- data.frame(data, p0, p1)

    # Normalize predictions:
    cols <- c('p0', 'p1')
    r_sum <- rowSums(data[cols])
    data$p0 <- data$p0 / r_sum
    data$p1 <- data$p1 / r_sum

    data <- data %>%
      mutate(pr_cx = p1,
             pr_ben = p0) %>%
      select(-c(p0, p1))
  }
  
  if (outcome == "ISUP"){
    # ISUP0
    lp0 <- with(data, log(X0 / (1 - X0)))
    fit0 <- glm(ISUP == 0 ~ lp0, family = binomial, data = data)
    p0 <- predict(fit0, type = "response")
    # ISUP1
    lp1 <- with(data, log(X1 / (1 - X1)))
    fit1 <- glm(ISUP == 1 ~ lp1, family = binomial, data = data)
    p1 <- predict(fit1, type = "response")
    # ISUP2
    lp2 <- with(data, log(X2 / (2 - X2)))
    fit2 <- glm(ISUP == 2 ~ lp2, family = binomial, data = data)
    p2 <- predict(fit2, type = "response")
    # ISUP3
    lp3 <- with(data, log(X3 / (3 - X3)))
    fit3 <- glm(ISUP == 3 ~ lp3, family = binomial, data = data)
    p3 <- predict(fit3, type = "response")
    # ISUP4
    lp4 <- with(data, log(X4 / (4 - X4)))
    fit4 <- glm(ISUP == 4 ~ lp4, family = binomial, data = data)
    p4 <- predict(fit4, type = "response")
    # ISUP5
    lp5 <- with(data, log(X5 / (5 - X5)))
    fit5 <- glm(ISUP == 5 ~ lp5, family = binomial, data = data)
    p5 <- predict(fit5, type = "response")
    
    # Updated Predictions with calibrated probabilities
    data <- data.frame(data, p0, p1, p2, p3, p4, p5)
    
    # Normalize predictions:
    cols <- c('p0','p1','p2','p3','p4','p5')
    r_sum <- rowSums(data[cols])
    data$p0 <- data$p0 / r_sum
    data$p1 <- data$p1 / r_sum
    data$p2 <- data$p2 / r_sum
    data$p3 <- data$p3 / r_sum
    data$p4 <- data$p4 / r_sum
    data$p5 <- data$p5 / r_sum
    
    data <- data %>%
      mutate(X0 = p0, X1 = p1, X2 = p2, X3 = p3, X4 = p4, X5 = p5) %>% 
      select(-c(p0:p5))
  }
  return(data)
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

