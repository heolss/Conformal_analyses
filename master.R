
rm(list=ls())

########################################################################################
## Main scripts
########################################################################################
# 1. Datasources - Table1
source("analysis/datasources.R")
datasources

# 2. Results - Baseline test set - Test-set 1
source("analysis/CP_LO_TEST.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 3. Results - Imagebase - Test set 2
source("analysis/CP_LO_IMAGEBASE.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 4. Results - External scanner and external pathology laboratory - Test set 4
source("analysis/CP_LO_KS.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 5. Results - CP on Rare prostate tissue morphology by subtype - Test set 6
source("analysis/CP_rare_subtypes_KS_by_subtype.R")
tab01;tab05;tab1
cx

# 6. Results - AI without CP on Rare prostate tissue morphology by subtype - Test set 6
source("analysis/AI_rare_subtypes_KS_by_subtypes.R")
AI_CX_subtype

# 8. Results - External scanner evaluted on 448 paired slides (Aperio & Hamamatsu) - Test set 3
source("analysis/paired_aperio.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)
source("analysis/paired_hamamatsu.R")
tab01;tab20;tab33
grid.arrange(cx, isup, ncol = 2)

# 9. Results - External scanner and external pathology laboratory - Test set 5
source("analysis/stavanger.R")
tab01;tab05;tab1
cx

# 10. Kolmogorov-Smirnov test of equality of the distribution of the predictions in the calibration set and 
# each test dataset to test the validity of the prediction regions. 
# The null hypothesis being that the samples are drawn from the same distribution.
source("analysis/kolmogorov_pvalues.R")
p_values

# 11. Prediction-regions by Imagebase-ISUP-mode and coverage by individual Imagebase-pathologist-votes
source("analysis/imagebase_voting_by_predgroups.R")
t1; t2; summary(df$Frequency)

################################## end of program #################################
