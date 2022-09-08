# Overview
Example project, providing demo versions of the analysis data and the scripts that are used for the analyses of the study “Estimating diagnostic uncertainty in artificial intelligence assisted pathology using conformal prediction”. Replicating the core tables and figures in the manuscript.

# Abstract
Unreliable predictions can occur when an artificial intelligence (AI) system is presented with data it has not been exposed to during training. We demonstrate the use of conformal prediction to detect unreliable predictions, using histopathological diagnosis and grading of prostate biopsies as example. We digitized 7788 prostate biopsies from 1192 men in the STHLM3 diagnostic study, used for training, and 1688 biopsies from 389 men used for testing. With conformal prediction, 1 in 794 (0.1%) predictions was incorrect while 175 (22%) of the predictions were flagged as unreliable when the AI-system was presented with new data from the same lab and scanner that it was trained on. Conformal prediction could with small samples (N=49 for external scanner and N=10 for external lab and scanner) detect systematic differences in external data. The AI-system with conformal prediction committed  3 (2%) errors for cancer detection in cases of atypical prostate tissue compared to 44 (25%) without conformal prediction, while flagging 143 (80%) unreliable predictions. We conclude that conformal prediction can increase patient safety of AI-systems.

# Repo-contents

## Data/Demo
We provide anonymized demo versions of all the datasets that are included in the manuscript. 

## Scripts

### functions.R
- Contains the main functions that are used for the main analysis of the paper

### master.R
- The analyses for the study are divided into different sub-scripts that are used to produce the different tables and figures of the manuscript.
- The master-script collects all sub-scripts together to provide an overview of the analysis. This cacan be used to run each sub-script / analysis in a serial way that aims to reproduce the structure of the analysis of how the results are presented in the manuscript.

## Output
- An output folder where tables and figures are collected. 

# System requirements
All statistical analyses were performed using R, version 4.0.1 (R Foundation). R-packages dplyr, tidyr, ggplo2, gridExtra and reshape2 were utilized for the analysis.
