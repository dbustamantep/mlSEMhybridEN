
# Multilevel ACES multivariate (based on bivariate model) for - TEs, sMRI volumne mediators, PTSDsxs 
# Using prepped data from H Maes and Daniel Bustamante: regressed out age, race, sex, scanner and whole brain volume

# ----------------------------------------------------------------------------------------------------------------------
# Program: mlevACES....R  
#  Author: Michael Neale, based on H Maes & M Hunter univariate; edited by Daniel Bustamante
#    Date: 11 27 2019 
#
# Purpose: FIT ACES model using MULTILEVEL Mediation modeling 
#
# -------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|

# Load Libraries & Options
rm(list=ls())

library(OpenMx)
library(psych)
library(matrixStats)
# library(car)
# library(umx) ##
source("miFunctions.R")
#source("miTryHard.R")
mxOption( NULL, "Default optimizer", "SLSQP" )
mxOption(NULL, 'Number of Threads', parallel::detectCores())
options(width=280)
options(max.print=999999)
mxVersion()


# Create Output 
filename    <- "mlMed_TEs_sMRIvol_PTSD"
sink(paste(filename,".Ro",sep=""), append=FALSE, split=TRUE)

# load data
nda18s <- readRDS("teMRIptsd_WholeBrain_Subcortical_data.Rds") 
# nda18c <- readRDS("teMRIptsd_WholeBrain_Cortical_data.Rds")

# select sMRI volume variables on indirect path
medVars <- grep(c("subcorticalgrayvolume|cerebellum_cortex_lh|cerebral_white_matter_lh|caudate_rh|caudate_lh|lateral_ventricle_rh"), medVars, value = TRUE)
medVars

predictor <- "rte"
outcome <- "rptsd_noTEna"
selVars <- c(  "rte", medVars, "rptsd_noTEna" )

# get variables with common names
# commonNames_nda18 <- intersect(names(nda18s), names(nda18c))
# bring the cortical variables into the subcortical dataset by IID
if (length(medVars_cort) > 0) {
nda18 <- merge(nda18s, nda18c[,c("IID", medVars_cort)], by = "IID")
} else {
  nda18 <- nda18s
}

rmuse <- selVars
rmuse

# Decide on whether to keep models and whether to scatterplot them - not recommended for large numbers of models
keepModels <- TRUE
splot <-  FALSE
# Decide whether to run everything against everything (slowest) or just the first one vs. rest 
firstVSrest  <- FALSE
malesOnly <- FALSE
femalesOnly <- FALSE # NB Making both malesOnly and femalesOnly TRUE will result in nobody at all. Not sure if 1 is Male yet(!!)
if(malesOnly & femalesOnly){stop("Not enough intersex people for analysis; set malesOnly or femalesOnly, not both")}
### End of change section: Rest shouldn't need changing


# ----------------------------------------------------------------------------------------------------------------------
numVars   <- c("relp","srelp","srelq","site","age","sex","race")
ndaz      <- nda18[,c("FID","IID",numVars,rmuse)] 
ndaz$FID <- as.factor(ndaz$FID)
ndaz$IID <- as.factor(ndaz$IID)
# Note, for MLM the join keys must be type factor

if(malesOnly){ndaz <- ndaz[ndaz$sex==2,]}
if(femalesOnly){ndaz <- ndaz[ndaz$sex==1,]}
# print a few summaries
dim(ndaz); head(ndaz); 
describe(ndaz[,numVars])
str(ndaz[,c("FID",numVars)])
ndazs     <- ndaz[order(ndaz$FID),]
table(ndaz$relp,  useNA="always")
table(ndaz$srelp, useNA="always")
table(ndaz$srelq, useNA="always")
table(ndaz$site,  useNA="always")



# ----------------------------------------------------------------------------------------------------------------------
# SHAPE DATA FOR MULTILEVEL MODELING

# Read tall data
tiDat       <- ndaz
dim(tiDat); names(tiDat)

# Delete individuals with missing definition variable ('relpi') for twins with unknown zygosity
tiDat       <- tiDat[!is.na(tiDat$relp),]
dim(tiDat); head(tiDat)

# ----------------------------------------------------------------------------------------------------------------------
# Make a nerdy amusing utility function or two to label matrices
labFun <- function(name="matrix",nrow=1,ncol=1,lower=FALSE,symmetric=FALSE,vector=F)
            {
                matlab <- matrix(paste(rep(name, each=nrow*ncol), rep(rep(1:nrow),ncol), rep(1:ncol,each=nrow),sep="_"),nrow,ncol)
				if((nrow > 1) | (ncol > 1))
								{
					                if(symmetric){for (i in 1:(nrow-1)){for (j in (i+1):ncol){matlab[i,j]<-matlab[j,i]}}}
                					if (lower) {for (i in 1:(nrow-1)){for (j in (i+1):ncol){matlab[i,j]<-NA}}}
									if(vector)
									{
										if(lower)
										{matlab<-as.vector(matlab[lower.tri(matlab)])} else{matlab<-as.vector(matlab)}
									}
								}
				return(matlab)
            }
lowerify  <- function(amatrix)
	{
		as.vector( amatrix[lower.tri((amatrix), diag=T)] )
	}
diagify  <- function(amatrix)
	{
		as.vector( amatrix[diag((amatrix))] )
	}

# ----------------------------------------------------------------------------------------------------------------------
# DEFINE MULTILEVEL MODELING FUNCTION

# Declare arguments for testing purposes

# Define function
runACESvcL <- function(dataset, varNames, famName, relpName, srelpName, srelqName, siteName, sexName, estimateCI=T, printStuff=T, tryHard=0, doc="none") 
{
if(printStuff){print(varNames)}
	fitError <- F; fitWarn <- F
	erFun 		<- function(e){print(paste("Error with variables",varNames));   fitError<-T}
	waFun 		<- function(w){print(paste("Warning with variables",varNames)); fitWarn <-T}
	finallyFun  <- function(f){if(fitWarn){fitMLM<-modelMLM}}

# ------------------------
# PREPARE MULTILEVEL DATA

# Select Dataset for Analysis
tallData  <- dataset[,c(varNames, famName, relpName, srelpName, srelqName, siteName, sexName)]
dim(tallData)

# Select Variables for Analysis
vars       <- varNames                   # list of variables' names (TODO check if it's a list)
nv         <- length(varNames)                         # number of variables
# It's rather annoying that colMeans/colVars doesn't work if there's only one column. Doh!
if(nv >  1){varMeans <- colMeans(tallData[,vars],na.rm=TRUE)}else{varMeans   <- mean(tallData[,vars],na.rm=TRUE)}
if(nv >  1){varVars  <- colVars (as.matrix(tallData[,vars]),na.rm=TRUE)}else{varVars   <- var(tallData[,vars],na.rm=TRUE)}

# Rescale variables (all) 
for(i in 1:nv){tallData[,vars[i]] <- (tallData[,vars[i]]-varMeans[i])/sqrt(varVars[i])}
# Select Data for Analysis:"within", "between" & "other" data sets
wData     <- tallData
bData     <- tallData[!duplicated(tallData[,famName]),  c(siteName, famName, relpName, srelpName, srelqName)]
oData     <- tallData[!duplicated(tallData[,siteName]), c(siteName, famName)]

if(printStuff){print(varMeans); print(varVars)}

# Set Starting Values
svMe      <- 0                   # start values for means
svPa      <- .2                # start values for variance component for A and C
svPc 	  <- svPa
svPe      <- .5		            # start values for variance component for E (TODO: check dimensions here)
svPs      <- .1                # start values for variance component for S
laCovAb   <- matrix(NA, nv, nv) 		# labels for diagonal of covAb
for (i in 1:nv)	{for (j in 1:nv) {laCovAb[i,j] <- paste0("relVa[",i,",",j,"]")}}
# -------------------------
# PREPARE MULTILEVEL MODEL

# Create Algebra for expected Mean Matrices
pathM     <- mxPath(from="one", to=vars, arrows=1, free=TRUE, values=svMe, labels=labVars("mean_",vars) )

# Create Unidirectional Paths
pathEw    <- mxPath( from=paste0('Ew',1:nv), to=vars, free=FALSE, values=1 )
pathAw    <- mxPath( from=paste0('Aw',1:nv), to=vars, free=FALSE, values=1, labels=paste0('data.',srelqName) )
pathAb    <- mxPath( from=paste0('between.','Ab',1:nv), connect='single', to=vars, free=FALSE, values=1, joinKey=famName ) 
pathCb    <- mxPath( from=paste0('between.','Cb',1:nv), connect='single', to=vars, free=FALSE, values=1, joinKey=famName )
pathOb    <- mxPath( from=paste0('between.','Ob',1:nv), connect='single', to=vars, free=FALSE, values=1, joinKey=famName )
pathOs    <- mxPath( from=paste0('other.'  ,'Os',1:nv), connect='single', to=paste0('Ob',1:nv), free=FALSE, values=rep(1,nv), joinKey=siteName )
#pathOs    <- mxPath( from=paste0('other.'  ,'Os',1:nv), connect='single', to=paste0('Ob',1:nv), free=FALSE, values=1, joinKey=siteName )
matrixP   <- mxMatrix( nrow=1, ncol=1, free=FALSE, values=1, labels=paste0('data.',relpName), name="P" )
matrixR   <- mxAlgebra( expression=P*Va, name="relVa" )

# Create CoVariance Paths
if(doc == "None")
	{
		covEw     <- mxPath( from=paste0('Ew',1:nv), to=paste0('Ew',1:nv), arrows=2, connect='unique.pairs', free=TRUE,  values=lowerify(diag(svPe,nv)), labels=lowerify(labFun('Ve',nv,nv,symmetric=T))    )
		covAw     <- mxPath( from=paste0('Aw',1:nv), to=paste0('Aw',1:nv), arrows=2, connect='unique.pairs', free=TRUE,  values=lowerify(diag(svPa,nv)), labels=lowerify(labFun('Va',nv,nv,symmetric=T))    )
		covAb     <- mxPath( from=paste0('Ab',1:nv), to=paste0('Ab',1:nv), arrows=2, connect='unique.pairs', free=FALSE, values=lowerify(diag(svPa,nv)), labels=lowerify(laCovAb) )
		covCb     <- mxPath( from=paste0('Cb',1:nv), to=paste0('Cb',1:nv), arrows=2, connect='unique.pairs', free=TRUE,  values=lowerify(diag(svPc,nv)), labels=lowerify(labFun('Vc',nv,nv,symmetric=T))    )
		covOb     <- mxPath( from=paste0('Ob',1:nv), to=paste0('Ob',1:nv), arrows=2, connect='unique.pairs', free=FALSE, values=0)
		covOs     <- mxPath( from=paste0('Os',1:nv), to=paste0('Os',1:nv), arrows=2, connect='unique.pairs', free=TRUE,  values=lowerify(diag(svPs,nv)), labels=lowerify(labFun('Vs',nv,nv,symmetric=T)) )
	}else if(doc == "AtoB" || doc == "AtoBnoBtwMeds") ## tell the function to make the cov paths if either AtoB or AtoBnoBtwMeds.
	{
		covEw     <- mxPath( from=paste0('Ew',1:nv), to=paste0('Ew',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPe,nv)), labels=lowerify(labFun('Ve',nv,nv,symmetric=T))    )
		covAw     <- mxPath( from=paste0('Aw',1:nv), to=paste0('Aw',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPa,nv)), labels=lowerify(labFun('Va',nv,nv,symmetric=T))    )
		covAb     <- mxPath( from=paste0('Ab',1:nv), to=paste0('Ab',1:nv), arrows=2, connect='unique.pairs', free=FALSE, 						   values=lowerify(diag(svPa,nv)), labels=lowerify(laCovAb) )
		covCb     <- mxPath( from=paste0('Cb',1:nv), to=paste0('Cb',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPc,nv)), labels=lowerify(labFun('Vc',nv,nv,symmetric=T))    )
		covOb     <- mxPath( from=paste0('Ob',1:nv), to=paste0('Ob',1:nv), arrows=2, connect='unique.pairs', free=FALSE                         ,  values=0)
		covOs     <- mxPath( from=paste0('Os',1:nv), to=paste0('Os',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPs,nv)), labels=lowerify(labFun('Vs',nv,nv,symmetric=T)) )
	}else if(doc == "BtoA")                                                                                                                               
	{                                                                                                                                                    
		covEw     <- mxPath( from=paste0('Ew',1:nv), to=paste0('Ew',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPe,nv)), labels=lowerify(labFun('Ve',nv,nv,symmetric=T))    )
		covAw     <- mxPath( from=paste0('Aw',1:nv), to=paste0('Aw',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPa,nv)), labels=lowerify(labFun('Va',nv,nv,symmetric=T))    )
		covAb     <- mxPath( from=paste0('Ab',1:nv), to=paste0('Ab',1:nv), arrows=2, connect='unique.pairs', free=FALSE, 						   values=lowerify(diag(svPa,nv)), labels=lowerify(laCovAb) )
		covCb     <- mxPath( from=paste0('Cb',1:nv), to=paste0('Cb',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPc,nv)), labels=lowerify(labFun('Vc',nv,nv,symmetric=T))    )
		covOb     <- mxPath( from=paste0('Ob',1:nv), to=paste0('Ob',1:nv), arrows=2, connect='unique.pairs', free=FALSE                         ,  values=0)
		covOs     <- mxPath( from=paste0('Os',1:nv), to=paste0('Os',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPs,nv)), labels=lowerify(labFun('Vs',nv,nv,symmetric=T)) )
	}else if(doc == "Both")                                                                                                                              
	{                                                                                                                                                    
		covEw     <- mxPath( from=paste0('Ew',1:nv), to=paste0('Ew',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPe,nv)), labels=lowerify(labFun('Ve',nv,nv,symmetric=T))    )
		covAw     <- mxPath( from=paste0('Aw',1:nv), to=paste0('Aw',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPa,nv)), labels=lowerify(labFun('Va',nv,nv,symmetric=T))    )
		covAb     <- mxPath( from=paste0('Ab',1:nv), to=paste0('Ab',1:nv), arrows=2, connect='unique.pairs', free=FALSE, 						   values=lowerify(diag(svPa,nv)), labels=lowerify(laCovAb) )
		covCb     <- mxPath( from=paste0('Cb',1:nv), to=paste0('Cb',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPc,nv)), labels=lowerify(labFun('Vc',nv,nv,symmetric=T))    )
		covOb     <- mxPath( from=paste0('Ob',1:nv), to=paste0('Ob',1:nv), arrows=2, connect='unique.pairs', free=FALSE                         ,  values=0)
		covOs     <- mxPath( from=paste0('Os',1:nv), to=paste0('Os',1:nv), arrows=2, connect='unique.pairs', free=as.logical(lowerify(diag(nv))),  values=lowerify(diag(svPs,nv)), labels=lowerify(labFun('Vs',nv,nv,symmetric=T)) )
	}else {stop(c("doc must be a character string AtoB, BtoA, Both, or None but you have it set to:",doc))}
	
## Causal paths for DOC Model  ##
Va <- mxMatrix("Full", nv, nv, free=as.logical(diag(nv)), labels=(labFun('Va',nv,nv,symmetric=T)), values=(diag(svPa,nv)), name="Va")
Vc <- mxMatrix("Full", nv, nv, free=as.logical(diag(nv)), labels=(labFun('Vc',nv,nv,symmetric=T)), values=(diag(svPc,nv)), name="Vc")
Ve <- mxMatrix("Full", nv, nv, free=as.logical(diag(nv)), labels=(labFun('Ve',nv,nv,symmetric=T)), values=(diag(svPe,nv)), name="Ve")
Vs <- mxMatrix("Full", nv, nv, free=as.logical(diag(nv)), labels=(labFun('Vs',nv,nv,symmetric=T)), values=(diag(svPs,nv)), name="Vs")
if(doc == 'AtoB')
	{
		phePath <-  mxPath( from=vars[1:(nv-1)],  to=vars[2:nv], connect='unique.bivariate', free=TRUE, values=0)
	}else
  if(doc == 'AtoBnoBtwMeds') ## doc with no causal paths between medVars
	{
	    from <- c(vars[1:(nv-1)], rep(vars[1], (nv-2))) ## get non-outcome  (pred and meds) variables for 'from' arg in mxPath (from 1st to all except outcome, plus as many predictors as mediators to be their A in the path; or pred, meds, rep pred # of meds)
	    to <- c(rep(vars[nv], nv-1), vars[2:(nv-1)])  ## get non-predictor (in this order: outcome and then med) variables for 'to' arg in mxPath (as many outcomes as the number of predictor and meds, then add the range of mediatoros (from the second var to the penultimate), plus the med vars as Bs of the predictor)
	    phePath <-  mxPath( from=from,  to=to, connect='single', free=TRUE, values=0) ## use 'single' in connect arg, to match their other with one AtoB, rather than all As to all Bs.
  }else
	if(doc == 'BtoA')
	{
		phePath <-  mxPath( from=vars[nv:2],  to=vars[(nv-1):1], connect='unique.bivariate', free=TRUE, values=0)
	}else
	if(doc == 'Both')
	{
		phePath <-  mxPath( from=vars, connect='all.bivariate', to=vars, free=TRUE, values=0) 
	}else
	if(doc == 'None')
	  
	{
		phePath <-  mxPath( from=vars, connect='unique.bivariate', to=vars, free=FALSE, values=0)
		Va <- mxMatrix("Full", nv, nv, free=T, labels=(labFun('Va',nv,nv,symmetric=T)), values=(diag(svPa,nv)), name="Va")
		Vc <- mxMatrix("Full", nv, nv, free=T, labels=(labFun('Vc',nv,nv,symmetric=T)), values=(diag(svPc,nv)), name="Vc")
		Ve <- mxMatrix("Full", nv, nv, free=T, labels=(labFun('Ve',nv,nv,symmetric=T)), values=(diag(svPe,nv)), name="Ve")
		Vs <- mxMatrix("Full", nv, nv, free=T, labels=(labFun('Vs',nv,nv,symmetric=T)), values=(diag(svPs,nv)), name="Vs")
		
		# Set path to nothing
	}else
	{stop(c("doc must be a character string AtoB, AtoBnoBtwMeds, BtoA, Both, or None but you have it set to:",doc))}


matrixV   <- mxAlgebra( expression=Va+Vc+Ve+Vs, name="V" )

# Create Data Objects for Multiple Groups
dataW     <- mxData( type="raw", observed=wData, sort=FALSE )
dataB     <- mxData( type="raw", observed=bData, primaryKey=famName )
dataO     <- mxData( type="raw", observed=oData, primaryKey= siteName )

# Create Model Objects for Multiple Groups
modelO    <- mxModel( 'other', type="RAM", dataO, latentVars=paste0('Os',1:nv), covOs )
modelB    <- mxModel( 'between', type="RAM", modelO, dataB, latentVars=c(paste0('Ab',1:nv),paste0('Cb',1:nv),paste0('Ob',1:nv)),
                       pathOs, covAb, covCb, covOb, matrixP, matrixR,  Va )
modelW    <- mxModel( 'within', type="RAM", modelB, dataW, manifestVars=vars, latentVars=c(paste0('Ew',1:nv),paste0('Aw',1:nv)),
                       pathEw, pathAw, pathAb, pathCb, pathOb, covEw, covAw, pathM, matrixV, Va, Vc, Ve, Vs )

# Create Algebra for Variance Components
rowV      <- paste0('estV',1:nv)
colV      <- paste0(rep(c('Va','Vc','Ve','Vs','Sa','Sc','Se','Ss'),each=nv),1:nv)
estV      <- mxAlgebra( expression=cbind(Va,Vc,Ve,Vs,Va/V,Vc/V,Ve/V,Vs/V), name="estV", dimnames=list(rowV,colV) )

# Create Confidence Interval Objects
if(doc=="None")
{
	ciACE     <- mxCI(paste0("estV[",eval(parse(text="1:nv")),",(4*nv+1):(8*nv)]" ))
	ciAlg <- NULL
}else{
	ciAlg <- mxAlgebra(cbind(A[2,1]*A[4,2],A[3,1]*A[4,3]),name="medEffects") ## make the algebra for medEffects paths adaptable for N medVars
	ciACE <- mxCI("medEffects")
}

# Build Model with Confidence Intervals
modelMLM <- mxModel( modelW, name='oneACEvcMLM', estV, ciACE, ciAlg, phePath ) 



# ---------------------
# RUN MULTILEVEL MODEL
fitMLM <- modelMLM

# Run ACE Model
## uncomment below (before commented out fitMLM)
# if(estimateCI)
# {
#   if(tryHard>0)
#   {fitMLM    <- mxPenaltySearchExternal(modelMLM)#mxTryHard( modelMLM, intervals=TRUE,  extraTries=tryHard, greenOK=T, silent=T)}else{ ## changed from mxTryHard to mxPenaltySearchExternal, mxregsem
#     (fitMLM <- mxPenaltySearchExternal(modelMLM)#mxRun( modelMLM, intervals=TRUE)
#       )} ## changed it from mxRun to mxPenaltySearchExternal, mxregsem
# } else{                                                                                                                                                          
#   if(tryHard<=0)                                                                                                                                               
#   {tryCatch(fitMLM    <- mxPenaltySearchExternal(modelMLM)#mxTryHard( modelMLM, intervals=FALSE, extraTries=0, greenOK=T, silent=T)
#             , error=erFun,  finally=finallyFun)}else{ ## changed from mxTryHard to mxPenaltySearchExternal, mxregsem
#     tryCatch(fitMLM <- mxPenaltySearchExternal(modelMLM)#mxRun( modelMLM, intervals=FALSE)
#              , error=erFun,  finally=finallyFun)} ## changed it from mxRun to mxPenaltySearchExternal, mxregsem
# }
if(estimateCI)
{
  if(tryHard>0)
  {fitMLM    <- mxTryHard( modelMLM, intervals=TRUE,  extraTries=tryHard, greenOK=T, silent=T)}else{ 
    (fitMLM <- mxRun( modelMLM, intervals=TRUE))}
} else{                                                                                                                                                          
  if(tryHard<=0)                                                                                                                                               
  {tryCatch(fitMLM    <- mxTryHard( modelMLM, intervals=FALSE, extraTries=0, greenOK=T, silent=T), error=erFun,  finally=finallyFun)}else{ 
    tryCatch(fitMLM <- mxRun( modelMLM, intervals=FALSE), error=erFun,  finally=finallyFun)}
}
#fitMLM    <- miTryHard(fitMLM, extraTries=2, greenOK=TRUE, bestInitsOutput=FALSE)

# Print Goodness-of-fit Statistics, Parameter Estimates & Confidence Intervals for Variance Components, if desired
if(printStuff)
{
fitGofs(fitMLM)
fitEsts(fitMLM)
if(estimateCI)
	{
		print((c(fitMLM$output$confidenceIntervals[1,1:nv],fitMLM$output$confidenceIntervals[2,1:nv],
                 fitMLM$output$confidenceIntervals[3,1:nv],fitMLM$output$confidenceIntervals[4,1:nv])))
	}
}              
# Create line with Stats & Estimates to be saved to output file              

npars <- length(fitMLM$output$estimate)
modsum <- summary(fitMLM)
if(estimateCI)
	{
		el   <-   c(round(modsum$de,0), fitMLM$output$min, (fitMLM$output$estimate), (fitMLM$output$standardErrors),  (fitMLM$output$confidenceIntervals), fitMLM$output$status$code, (modsum$cpu)  ,varNames)
#						el   <-   c(summary(fitMLM)$de, fitMLM$output$min, (fitMLM$output$estimate),  (fitMLM$output$standardErrors), fitMLM$output$confidenceIntervals, fitMLM$output$status$code, fitMLM$output$cpuTime, varNames)
  		 	print(names(el))
# TODO: fix names on confidence intervals!
			names(el)[1:2] <- c("df","neg2LL")
#  		 	dimnames(el)[(4+2*npars):(4+3*npars-1)] <- paste0("LowCI_",names(fitMLM$output$estimate))# needs work to get low, est, high order etc
			names(el)[(length(el)-nv+1):(length(el))] <- paste0("Var",1:nv)
	}else{if(length(fitMLM$output)>0)
			{
		
			el   <-   c(round(modsum$de,0), fitMLM$output$min, (fitMLM$output$estimate), (fitMLM$output$standardErrors), fitMLM$output$status$code, (modsum$cpu), varNames)
#			el   <-   c(summary(fitMLM)$de, fitMLM$output$min, (fitMLM$output$estimate),  (fitMLM$output$standardErrors), fitMLM$output$status$code, fitMLM$output$cpuTime, varNames)
			}else{el <- c(rep(NA,2+2*(4*nv*(nv+1)/2)),varNames)}
   		 	print(names(el))
   		 	names(el)[1:2] <- c("df","neg2LL")
   		 	names(el)[(3+npars):(3+2*npars-1)] <- paste0("SE_",names(fitMLM$output$estimate))
			names(el)[(3+2*npars):(3+2*npars+1)] <- c("status","cpu")
			names(el)[(3+2*npars+2):(3+2*npars+2+nv-1)] <- paste0("Var",1:nv)
		 }
	 
return(list(el=el,fit=fitMLM))

# End of function
}	


abcdVars     <- c(rmuse)
abcdLength  <- length(abcdVars)

docString <- "AtoBnoBtwMeds"

if(keepModels){xMLM        <- vector(mode="list", length=abcdLength*(abcdLength-1)/2)}

### Just run one multivariate model with all of rmuse (i.e. abcdVars) list
justOne <- T
if(justOne){
ij <- 1
  if(splot)
{
  {pdf("noPihat.pdf")}
  mystr <- paste0("~",abcdVars[i],"+",abcdVars[j])
  scatterplotMatrix(as.formula(mystr), data=tiDat)
}
if(keepModels)
{
  xMLM[[ij]] <- runACESvcL(dataset=tiDat, varNames=c(abcdVars), famName='FID', relpName='relp', srelpName='srelp', srelqName='srelq', siteName='site', sexName='sex', estimateCI=F, printStuff=F, tryHard=2, doc=docString)
  print(xMLM[[ij]]$el)
  if (ij == 1) { out <- xMLM[[ij]]$el }else{ out <- rbind(out, unname(xMLM[[ij]]$el))}
}else
{
  x <- runACESvcL(dataset=tiDat, varNames=c(abcdVars), famName='FID', relpName='relp', srelpName='srelp', srelqName='srelq', siteName='site', sexName='sex', estimateCI=F, printStuff=F, tryHard=2, doc=docString)
  print(x$el)
  if (ij == 1) { out <- x$el }else{ out <- rbind(out, unname(x$el))}
}
write.table(out,file=paste(filename,".csv",sep=""),sep=",")
}

# save first fit:
modelFit <- xMLM[[1]]$fit
save(modelFit, file = paste(filename, "N", nrow(nda18), length(medVars), "meds_1stFit.RData", sep = "_")) 
summary(xMLM[[1]]$fit)
summary(modelFit)

modelFit$oneACEvcMLM

# Modify model to allow bivariate ACE for mediatiors, but causal otherwise
xMLM[[1]]$fit$A$free
xMLM[[1]]$fit$S

modelFit2 <- omxSetParameters(modelFit, labels=c("Va_5_4","Vc_5_4","Ve_5_4", "Vs_5_4"
                                                 ), free=T) 
modelFit2 <- mxModel(modelFit2, name = "modelFit2_ROI_LR_aces") 
modelfit2out <- mxTryHard(modelFit2, intervals = FALSE)
summary(modelfit2out)

# uncomment to run model comparisons
# # TEs <-> PTSDsxs
# biDirPath <- modelfit2out
# biDirPath <- mxModel(biDirPath, name = "biDirCpath")
# biDirPath$A$free[1,(length(rmuse))] <- TRUE
# biDirPathOut <- mxTryHard(biDirPath, intervals = FALSE)
# summary(biDirPathOut)
# 
# comp1 <- mxCompare(biDirPathOut, list(modelfit2out, modelFit))
# comp1
# 
# save(biDirPathOut, file = paste(filename, "N", nrow(nda18), length(medVars), "meds_biCpathFit.RData", sep = "_")) 
#  

# # TEs <-> PTSDsxs and MedL-R <-> PTSDsxs
# biPaths2 <- biPathsout
# biPaths2$A$free[2,4] <- TRUE
# biPaths2$A$free[3,4] <- TRUE
# biPath2sout <- mxTryHard(biPaths2, intervals = TRUE)
# summary(biPath2sout) 


# # TEs <-> PTSDsxs and TEs <-> MedL-R <-> PTSDsxs
# allBiPaths <- biDirPathOut
# allBiPaths <- mxModel(allBiPaths, name = "allBiPaths")
# allBiPaths$A$free[1,(2:(length(medVars)+1))] <- TRUE
# allBiPaths$A$free[(2:(length(medVars)+1)), length(rmuse)] <- TRUE
# allBiPathsOut <- mxTryHard(allBiPaths, intervals = FALSE)
# summary(allBiPathsOut)
# 
# comp2 <- mxCompare(allBiPathsOut, list(biDirPathOut, modelfit2out, modelFit))
# comp2
# 
# 
# allBiPathsOut$A
# allBiPathsOut$S
# 
# save(allBiPathsOut, file = paste(filename, "N", nrow(nda18), length(medVars), "meds_allBiPathsFit.RData", sep = "_")) ## gsub('.{3}$', '', rmuse[2]) for gsub pattern argument: "." any character, "{n}" matched n times, "$"  at the end of line




# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
sink()
#save.image(paste(filename,".Ri",sep=""))
save.image(file = paste(filename,"N", nrow(nda18), length(medVars),"meds_workspace.RData", sep = "_"))

