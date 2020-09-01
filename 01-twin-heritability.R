# -------------------------------------------------------------------------
# Estimate heritablity for microbiome taxa from twin data
# -------------------------------------------------------------------------




# SETUP -------------------------------------------------------------------
library(OpenMx)
source("R/00_miFunctions.R")      # Load helper functions
mxVersion()




# LOAD DATA ---------------------------------------------------------------
# Taxa proportion data should be in separate objects for MZ and DZ pairs with
# 2 columns and N rows corresponding to the number of pairs.
# Column 1 contains taxa data for the first twin and column 2 contains data
# for the second twin of the same pair.
load("path/to/data")





# SET MODEL PARAMETERS ----------------------------------------------------

vars      <- 'taxa'                                      # list of variable names
nv        <- 1                                           # number of variables
ntv       <- nv*2                                        # number of total variables
selVars   <- paste(vars,c(rep(1,nv),rep(2,nv)),sep="")


# Adjust Starting Values
svMe      <- .25                                         # start value for means
svPa      <- .2                                          # start value for path coefficient
svPe      <- .5                                          # start value for path coefficient for e


# Create Algebra for expected Mean Matrices
meanG     <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean",vars), name="meanG" )

# Create Matrices for Variance Components
covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VA11", name="VA" )
covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VC11", name="VC" )
covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=svPa, label="VE11", name="VE" )

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP      <- mxAlgebra( expression= VA+VC+VE, name="V" )
covMZ     <- mxAlgebra( expression= VA+VC, name="cMZ" )
covDZ     <- mxAlgebra( expression= 0.5%x%VA+ VC, name="cDZ" )
expCovMZ  <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ  <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )

# Create Data Objects for Multiple Groups
dataMZ    <- mxData( observed=mzData, type="raw" )
dataDZ    <- mxData( observed=dzData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ     <- mxExpectationNormal( covariance="expCovMZ", means="meanG", dimnames=selVars )
expDZ     <- mxExpectationNormal( covariance="expCovDZ", means="meanG", dimnames=selVars )
funML     <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
pars      <- list( meanG, covA, covC, covE, covP )
modelMZ   <- mxModel( pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ   <- mxModel( pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
multi     <- mxFitFunctionMultigroup( c("MZ","DZ") )

# Create Algebra for Unstandardized and Standardized Variance Components
rowUS     <- rep('US',nv)
colUS     <- rep(c('VA','VC','VE','SA','SC','SE'),each=nv)
estUS     <- mxAlgebra( expression=cbind(VA,VC,VE,VA/V,VC/V,VE/V), name="US", dimnames=list(rowUS,colUS) )

# Create Confidence Interval Objects
ciACE     <- mxCI( "US[1,1:3]" )

# Build Model with Confidence Intervals
modelACE  <- mxModel( "oneACEvc", pars, modelMZ, modelDZ, multi, estUS, ciACE )


# Run ACE Model
fitACE    <- mxRun( modelACE, intervals = TRUE, useOptimizer = TRUE )
sumACE    <- summary( fitACE )

# Compare with Saturated Model
#if saturated model fitted in same session
#mxCompare( fitSAT, fitACE )
#if saturated model prior to genetic model
#lrtSAT(fitACE,4055.93‚46,1767)

# Print Goodness-of-fit Statistics & Parameter Estimates
# fitGofs(fitACE)
# fitEstCis(fitACE)

# Run AE model
modelAE   <- mxModel( fitACE, name="oneAEvc" )
modelAE   <- omxSetParameters( modelAE, labels="VC11", free=FALSE, values=0 )
fitAE     <- mxRun( modelAE, intervals=T , useOptimizer = TRUE )
# fitGofs(fitAE); fitEstCis(fitAE)

# Run CE model
modelCE   <- mxModel( fitACE, name="oneCEvc" )
modelCE   <- omxSetParameters( modelCE, labels="VA11", free=FALSE, values=0 )
modelCE   <- omxSetParameters( modelCE, labels=c("VE11","VC11"), free=TRUE, values=.6 )
fitCE     <- mxRun( modelCE, intervals=T ,useOptimizer = TRUE )
# fitGofs(fitCE); fitEstCis(fitCE)

# Run E model
modelE    <- mxModel( fitAE, name="oneEvc" )
modelE    <- omxSetParameters( modelE, labels="VA11", free=FALSE, values=0 )
fitE      <- mxRun( modelE, intervals=T ,useOptimizer = TRUE )
# fitGofs(fitE); fitEstCis(fitE)

# Print Comparative Fit Statistics

return(list(
  fitEstCis(fitAE),
  #fitEstCis(fitACE),
  mxCompare( fitACE, nested <- list(fitAE, fitCE, fitE) ),
  round(rbind(fitACE$US$result,fitAE$US$result,fitCE$US$result,fitE$US$result), 4)
))
