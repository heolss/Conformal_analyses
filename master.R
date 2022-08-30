
rm(list=ls())

########################################################################################
## Main scripts
########################################################################################
# 1. Datasources - Table1
source("datasources.R")
datasources

# 2. Results - Test-set
source("CP_LO_argmaxvoting/CP_LO_TEST.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 3. Results - Imagebase
source("CP_LO_argmaxvoting/CP_LO_IMAGEBASE.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 4. Results - KS
source("CP_LO_argmaxvoting/CP_LO_KS.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 5. Results - CP on Rare prostate tissue morphology by subtype
source("rare_subtypes/CP_rare_subtypes_KS_by_subtype.R")
tab01;tab05;tab1
cx

# 6. Results - AI without CP on Rare prostate tissue morphology by subtype
source("rare_subtypes/AI_rare_subtypes_KS_by_subtypes.R")
AI_CX_subtype

# 8. Results - STHLM3-Aperio model evaluted on 448 paired slides (scanned on both Aperio & Hamamatsu)
source("STHLM3_APERIO/paired_aperio.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)
source("STHLM3_APERIO/paired_hamamatsu.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 9. Prediction-regions by Imagebase-ISUP-mode and coverage by individual Imagebase-pathologist-votes
source("Programs/CP_LO_argmaxvoting/imagebase_voting_by_predgroups.R")
t1; t2; summary(df$Frequency)

# 10. Kolmogorov-Smirnov test of equality of the distribution of the predictions in the calibration set and 
# each test dataset to test the validity of the prediction regions. 
# The null hypothesis was that the samples are drawn from the same distribution.
source("CP_LO_argmaxvoting/kolmogorov_pvalues.R")
p_values


################################## end of program #################################
