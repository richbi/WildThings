rm(list=ls())
library(nimble)
library(raster)

home<-getwd()

##########
## BEAR ##
##########
## LOAD NIMBLE FUNCTIONS AND DATA 



setwd(file.path(home,"ScriptAndDataMCMC/Bear"))


## LOAD FEMALE DATA
load("17.F_12_18_INPUTChain1.RData")
## LOAD MALE DATA 
#load("17.M_12_18_INPUTChain1.RData")

## LOAD THE CUSTOM NIMBLE FUNCTIONS
source("dbin_LESSCachedAllSparseBear_v2.R")
source("pointProcess.R")

#RUN NIMBLE MODEL (demonstration with single chain and low number of iterations)
model <- nimbleModel( code = modelCode
                      , constants = nimConstants
                      , data = nimData
                      , inits = nimInits
                      , check = FALSE       
                      , calculate = FALSE)  
cmodel <- compileNimble(model)
cmodel$calculate()
MCMCconf <- configureMCMC(model = model, monitors = c(nimParams),
                          control = list(reflective = TRUE, adaptScaleOnly = TRUE),
                          useConjugacy = FALSE) 
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
Runtime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC
                                                  , nburnin = 0
                                                  , niter = 100
                                                  , nchains = 1
                                                  , samplesAsCodaMCMC = TRUE))

#--PLOT STORED ANNUAL DENSITY RASTERS (average posterior utilization density, see Methods) 

setwd(file.path(home,"DensityRasterMaps"))
load("DensityRasterBrickBear.RData")
plot(DensityRasterBrick)



##########
## WOLF ##
##########


setwd(file.path(home,"ScriptAndDataMCMC/Wolf"))


## LOAD FEMALE DATA
load("9.F1218Cached_INPUTChain1.RData")
## LOAD MALE DATA 
#load("9.M1218Cached_INPUTChain1.RData")

## LOAD THE CUSTOM NIMBLE FUNCTIONS
source("dbin_LESSCachedAllSparseWolf.R")
source("pointProcess.R")

#RUN NIMBLE MODEL (demonstration with single chain and low number of iterations)
model <- nimbleModel( code = modelCode
                      , constants = nimConstants
                      , data = nimData
                      , inits = nimInits
                      , check = FALSE       
                      , calculate = FALSE)  
cmodel <- compileNimble(model)
cmodel$calculate()
MCMCconf <- configureMCMC(model = model, monitors = c(nimParams),
                          control = list(reflective = TRUE, adaptScaleOnly = TRUE),
                          useConjugacy = FALSE) 
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
Runtime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC
                                                  , nburnin = 0
                                                  , niter = 100
                                                  , nchains = 1
                                                  , samplesAsCodaMCMC = TRUE))

#--PLOT STORED ANNUAL DENSITY RASTERS (average posterior utilization density, see Methods) 

setwd(file.path(home,"DensityRasterMaps"))
load("DensityRasterBrickWolf.RData")
plot(DensityRasterBrick)




###############
## WOLVERINE ##
###############


setwd(file.path(home,"ScriptAndDataMCMC/Wolverine"))

## LOAD FEMALE DATA
load("22.J_Fa1.RData")
## LOAD MALE DATA 
#load("22.J_Ma1.RData")

## LOAD THE CUSTOM NIMBLE FUNCTIONS
source("dbin_LESS_Cached_MultipleCovResponse.R")
source("pointProcess.R")

#RUN NIMBLE MODEL (demonstration with single chain and low number of iterations)
model <- nimbleModel( code = modelCode
                      , constants = nimConstants
                      , data = nimData
                      , inits = nimInits
                      , check = FALSE       
                      , calculate = FALSE)  
cmodel <- compileNimble(model)
cmodel$calculate()
MCMCconf <- configureMCMC(model = model, monitors = c(nimParams),
                          control = list(reflective = TRUE, adaptScaleOnly = TRUE),
                          useConjugacy = FALSE) 
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
Runtime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC
                                                  , nburnin = 0
                                                  , niter = 100
                                                  , nchains = 1
                                                  , samplesAsCodaMCMC = TRUE))

#--PLOT STORED ANNUAL DENSITY RASTERS (average posterior utilization density, see Methods) 

setwd(file.path(home,"DensityRasterMaps"))
load("DensityRasterBrickWolverine.RData")
plot(DensityRasterBrick)