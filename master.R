########################################################################################
## Read data and load packages
########################################################################################
rm(list=ls())

########################################################################################
## Set parameters
########################################################################################

########################################################################################
## Main scripts
########################################################################################
# 1. Datasources - Table1
source("Programs/datasources.R")
datasources

# 2. Results - Test-set
source("Programs/CP_LO_argmaxvoting/CP_LO_TEST.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 3. Results - Imagebase
source("Programs/CP_LO_argmaxvoting/CP_LO_IMAGEBASE.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 4. Results - CP on Rare prostate tissue morphology (Lars-special-cases)
source("Programs/CP_LO_argmaxvoting/CP_LO_LARS.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 5. Results - CP on Rare prostate tissue morphology by subtype (Lars-special-cases)
source("Programs/CP_LO_argmaxvoting/CP_LO_LARS_SUBTYPE.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 7. Results - AI without CP on Rare prostate tissue morphology by subtype (Lars-special-cases)
source("Programs/CP_LO_argmaxvoting/AI_LARS_SUBTYPES.R")
AI_CX_subtype

# 8. Results - STHLM3-Aperio model evaluted on 448 paired slides (scanned on both Aperio & Hamamatsu)
source("Programs/STHLM3_APERIO/paired_aperio.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)
source("Programs/STHLM3_APERIO/paired_hamamatsu.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)
source("Programs/STHLM3_APERIO/paired_calibration.R")

# 6. Prediction-regions by Imagebase-ISUP-mode and coverage by individual Imagebase-pathologist-votes
source("Programs/CP_LO_argmaxvoting/imagebase_voting_by_predgroups.R")
t1; t2; summary(df$Frequency)

########################################################################################
## Extra
########################################################################################
# Prediction-regions by Imagebase-consensus
source("Programs/CP_LO_argmaxvoting/CP_LO_IMAGEBASE_CONSENSUS.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

########################################################################################
## Notes
########################################################################################