
# ----------------------------------------------------------------------------------------------------------------------
#  Author: Daniel Bustamante
#    Date: 02 24 2020 
#
# Purpose: fit mediation model using mxregsem for a direct and indirect effects model (predictor, variables in the indirect path, outcome). Full sample
#
# ----------------------------------------------------------------------------------------------------------------------

# Load Libraries & Options
library(OpenMx)
library(psych)
library(matrixStats)
library(MASS)

## uncomment these two lines below if mxregem has not been installed --------- #
# install.packages("devtools")
# devtools::install_github("trbrick/mxregsem")
# library(mxregsem)
## --------------------------------------------------------------------------- #

source("miFunctions.R")
#source("miTryHard.R")
mxOption( NULL, "Default optimizer", "SLSQP" )
mxOption(NULL, 'Number of Threads', parallel::detectCores())
options(width=280)
options(max.print=999999)
options(scipen = 999)


# Create Output 
filename    <- "reg_subc"

# Load data
nda201 <- readRDS("data.Rds")


dim(nda201)
abcdNames<-names(nda201)
abcdNames[1:10]

## if variables' names have special characters, the change them, for example "." in the variables' names, for a "_" to avoid error while running models 
colnames(nda201) <- gsub("\\.+","\\_", colnames(nda201)) ## find any punctuation or space character as pattern, "+" (the preeceding items can be matched) one or more times, with \\ to treat "_" as a regular character in the replacement, could also do "[[:punct:][:space:]]+" if punctuation signs are unknown

# check some correlations
someVars <- c("rsmri_vol_subcort_aseg_hippocampus_lh", "rsmri_vol_subcort_aseg_hippocampus_rh",
              "rsmri_vol_subcort_aseg_amygdala_lh", "rsmri_vol_subcort_aseg_amygdala_rh",
              "rsmri_vol_subcort_aseg_putamen_lh", "rsmri_vol_subcort_aseg_putamen_rh",
              "rsmri_vol_subcort_aseg_caudate_lh", "rsmri_vol_subcort_aseg_caudate_rh")

cor_vars_subc <- nda201[, c("rte", "rptsd_noTEna", 
                            someVars
)]

round(cor(cor_vars_subc, use = "complete.obs"), 4)
cor.test(nda201[,someVars[1]], nda201[,someVars[2]])

# describe/check future mediation variables
nda201sMRIdesc <- describe(nda201[, grep(c("rsmri_vol_subcort_aseg"), names(nda201), value = TRUE)])
nda201sMRIdesc

# if any variable is full of NAs (meeting any of the conditions below), then identify them
naMedVars <- which(nda201sMRIdesc$n == 0 | is.nan(nda201sMRIdesc$mean) | is.na(nda201sMRIdesc$sd) | nda201sMRIdesc$range == -Inf)
nopeMedVars <- rownames(nda201sMRIdesc[naMedVars[1:length(naMedVars)],])
nopeMedVars

## select ROIs to include as mediators to be regularized (subcort, cort, etc)
medVars <- grep(c("rsmri_vol_subcort_aseg"), names(nda201), value = TRUE) ## getting ROIs for analysis 

# select ROIs to include as mediators
# remove variables not passing the conditions above
if (length(nopeMedVars) > 0) {medVars <- medVars[!grepl(paste0(nopeMedVars, collapse = "|"), medVars)]} else {medVars} ## with the "|" for using all the names in the object
medVars

# select ROIs to include as mediators; and predictor and outcome
predictor <- "rte"
outcome <- "rptsd_noTEna"

# data into analysis
rmuse <- c(predictor, medVars, outcome)
 
describe(nda201[, rmuse])



# ---------------------------------------
# Prepare variables and starting values:
pred <- predictor
med <- medVars
outc <- outcome

vars <- rmuse
manVars <- c(pred, med, outc)
nVars <- length(manVars)
sVa <- .1

modelDirection <- "TEs_volSubc_ptsd"

# run regularization?
regularize <- TRUE

# select type of regularization (ridge, LASSO, elastic net)
LASSO <- FALSE
ElasticNet <- TRUE

# select position to start in sequence for the penalty value if running reg model
# penPos <- 11 ## as the penalty val increases (passed to lambda), stronger penalty, then greater shrinkage of coefficients
lambda_minEN <- 0.025181323 
lambda_maxEN <- 0.12345954 # the largest lambda that keeps MSE within one SE of the minimal MSE
lambda_minLASSO <- 0.029959404 
lambda_maxLASSO <- 0.043465966 # the largest lambda that keeps MSE within one SE of the minimal MSE

# alpha is for the elastic-net mixing parameter, with range 0:1, α=1 is the lasso (default) and α=0 is the ridge
alpha_EN <- .8


# ---------------------------------------
# Prepare objects and model:

# paths
directPath <- mxPath(from = pred, to = outc, arrows = 1, labels = "c_dir", free = TRUE, values = sVa)
preMedPath <- mxPath(from = pred, to = med, arrows = 1, labels = paste("a_path", 1:length(med), med, sep = "_"), free = TRUE, values = sVa)
postMedPath <- mxPath(from = med, to = outc, arrows = 1, labels = paste("b_path", 1:length(med), med, sep = "_"), free = TRUE, values = sVa)

paths <- list(directPath, preMedPath, postMedPath)

# variance
preVar <- mxPath(from = pred, arrows = 2, free = TRUE, values = 1, labels="predVar")
medVar <- mxPath(from = c(med), arrows = 2, free = TRUE, values = 1, labels = paste("medVar_", c(med), sep = ""))
outcVar <- mxPath(from = outc, arrows = 2, free = TRUE, values = 1, labels="outcVar")

varcovs <- list(preVar, medVar, outcVar)

# means, data type, fit function and CIs
means <- mxPath(from = "one", to = c(manVars), free = TRUE, labels = paste("mean_", c(manVars), sep = ""))
dataMed <- mxData(nda201, type="raw")
funML <- mxFitFunctionML()
# ci <- mxCI(c("med"))

# build model including CIs
modelMed <- mxModel(paste("med_model", modelDirection, sep = "_"), type = "RAM", manifestVars = manVars,
                    paths, varcovs, means, dataMed, funML)#, ci)


# ---------------------------------------
if (regularize == TRUE) {
  library(mxregsem) ## call library here so mxregsem mxModel() doesn't mask OpenMx mxModel() creating and creates a _Reg model with PenaltyFunction when Regularization is not used
  # create mxregsem regularization object before model ----------------------------------------------------- ########
  if (LASSO == TRUE){
    # pen.seq = c(0,(seq(1,sqrt(nrow(nda18)),length.out=99)**2)) 
    # pen.val <- pen.seq[penPos] ## as the penalty val increases (passed to lambda), stronger penalty, then greater shrinkage of coefficients -- ###
  tLasso <- mxRegularizeLASSO(what = getParamsInMatrix(modelMed, "A"), name = "LASSO", 
                              lambda = lambda_minLASSO, lambda.max = lambda_maxLASSO)#, lambda.step = 1)
  
  modelMed <- mxModel(modelMed, tLasso) ## add regularization object tLasso to mxModel mxregsem
  }
  if (ElasticNet == TRUE){
    tElasticN <- mxRegularizeElasticNet(what = getParamsInMatrix(modelMed, "A"), name = "ElasticNet",
                                        alpha = alpha_EN, 
                                lambda = lambda_minEN, lambda.max = lambda_maxEN)
    
    modelMed <- mxModel(modelMed, tElasticN) ## add regularization object tLasso to mxModel mxregsem
    pen <- 0
    pen.val <- 0
    lambda_minLASSO <- 0
    lambda_maxLASSO <- 0 
  }
# Run reg model
fitMed  <- mxPenaltySearchExternal(modelMed, intervals = FALSE)

# ------------------------------------------------------------------------------------------------------- ########
# run non-regularized model
} else {
  fitMed  <- mxTryHard(modelMed, intervals = FALSE)
  # set penalty stuff to 0 if did not run regmodel so output can be saved with info anyway
  pen <- 0
  pen.val <- 0
  lambda_minEN <- 0
  lambda_maxEN <- 0
  lambda_minLASSO <- 0
  lambda_maxLASSO <- 0 
  }
# end of non-regularized or regularized run area
summary(fitMed)


if (ElasticNet == TRUE && LASSO == FALSE){
saveRDS(fitMed, file = paste(filename, "Direction", modelDirection, "meds", length(med), "N", nrow(nda201), "reg", regularize, "ElasticNet", ElasticNet, "LASSO", LASSO, "LambdaMin", round(lambda_minEN,4), "LambdaMax", round(lambda_maxEN,4), "alpha", fitMed$output$regularizationSearchOutcomes$ElasticNet_alpha[1], ".Rds", sep = "_" ))
} else {
saveRDS(fitMed, file = paste(filename, "Direction", modelDirection, "meds", length(med), "N", nrow(nda201), "reg", regularize, "ElasticNet", ElasticNet, "LASSO", LASSO, "LambdaMin", round(lambda_minLASSO,4), "LambdaMax", round(lambda_maxLASSO,4), ".Rds", sep = "_" ))
}


# ----------------------------------------------------------------------------------------------------------------------
save.image(file = paste(filename, "meds", length(med), "N", nrow(nda201), "reg", regularize, "pen_val", round(pen.val,4), "reg_pen", round(pen,4), "meds_workspace.RData", sep = "_"))


