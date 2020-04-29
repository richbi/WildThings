## 0. ------ DEREGISTER ANY PREVIOUSLY REGISTERED DISTRIBUTIONS ------
# In some versions of NIMBLE the distribution registeration proceedures can return an error if the distributions
# have already been previously registered.  This little section of code deregisters the functions to ensure no
# errors are thrown if this file is sourced more than once.
if(exists("distributions", nimbleUserNamespace)) {
   # List of distributions defined in this source file
   distributionNames <- c(
      "dbinomPP",
      "dbinomPPSingle",
      "dbinomMNormSourcePP",
      "dbinomMNormSourcePPSingle",
      "dbinomMNormSourcePPMulti",
      "dpoisPP",
      "dpoisMNormSourcePP",
      "dbinomMNormSourcePP_dbinomPP_SingleMixture")
   # Find those distributions that are defined in this source fiel and see if they are already registered
   isDefinedDist <- distributionNames %in% nimbleUserNamespace$distributions$namesVector
   if(any(isDefinedDist)) {
      # If the distributions are registered then deregister them
      deregisterDistributions(distributionNames[isDefinedDist])
   }
}

## 1. ------ DEFINE NIMBLE CUSTOM UTILITY FUNCTIONS USED IN THE POINT PROCESS DISTRIBUTIONS ------

### 1.1. ==== Calcualte the size of the observation window ====
calcObsWindowSize <- nimbleFunction(
   run = function(
      lowerCoords = double(2),
      upperCoords = double(2),
      areAreas = double(0, default = 1)
   ) {
      ## 1.1.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(1))
      ## 1.1.2. Sanity test the inputs ----
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(lowerCoords)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      }
      # Assess the number of observation windows
      numObsWindows <- dim(lowerCoords)[1]
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      } else if(numObsWindows != dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates")
      }
      if(dimCoords != dim(upperCoords)[2]) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      ## 1.1.3. Calculate the area/length of the observation window ----
      obsWindowSize <- rep(0, numObsWindows)
      if(areAreas) {
         # The observation windows are areas/volumes and therefore calculate their area/volume
         # accordingly
         obsWindowSize <- obsWindowSize + 1.0
         for(dimIter in 1:dimCoords) {
            obsWindowSize <- obsWindowSize * (upperCoords[1:numObsWindows, dimIter] - lowerCoords[1:numObsWindows, dimIter])
            if(sum(obsWindowSize <= 0.0) > 0) {
               # Ensure that upper and lower coordinates are valid
               stop("window area is zero for at least one observation window")
            }
         }
      } else {
         # The observation windows are transects and therefore calculate their length accordingly
         for(dimIter in 1:dimCoords) {
            coordDist <- upperCoords[1:numObsWindows, dimIter] - lowerCoords[1:numObsWindows, dimIter]
            obsWindowSize <- obsWindowSize + coordDist * coordDist
         }
         if(sum(obsWindowSize <= 0.0) > 0) {
            # Ensure that upper and lower coordinates are valid
            stop("window length is zero for at least one observation window")
         }
         obsWindowSize <- sqrt(obsWindowSize)
      }
      return(obsWindowSize)
   }
)

### 1.2. ==== Calculate the normalisation integral of the multivariate normal source-point model ====
calcMNormSourceInt <- nimbleFunction(
   run = function(
      lowerCoords = double(2),
      upperCoords = double(2),
      sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      areAreas = double(0, default = 1)
   ) {
      ## 1.2.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(1))
      ## 1.2.2. Sanity test the inputs ----
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(lowerCoords)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      }
      # Assess the number of observation windows
      numObsWindows <- dim(lowerCoords)[1]
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      } else if(numObsWindows != dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates")
      }
      if(dimCoords != dim(upperCoords)[2]) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Ensure that the source coordinates have the correct dimensionality
      if(length(sourceCoords) != dimCoords) {
         stop("source coordinates do match the coordinates of the observation windows")
      }
      # Ensure that the standard deviation of the isotropic multivariate normal distribution
      # has a valid value
      if(normSD <= 0.0) {
         stop("invalid value given for the standard deviation of the isotropic normal distribution")
      }
      ## 1.2.3. Calculate the normalisation integral ----
      # Initialise a vector to store the output
      normIntegral <- numeric(length = numObsWindows, value = 0.0, recycle = TRUE)
      if(areAreas) {
         # The observation windows are areas/volumes
         normIntegral <- normIntegral + pow(2.0 * pi * normSD * normSD, dimCoords / 2.0)
         # Iterate over the number of dimensions
         for(dimIter in 1:dimCoords) {
            if(sum(upperCoords[1:numObsWindows, dimIter] <= lowerCoords[1:numObsWindows, dimIter]) > 0) {
               # Ensure that upper and lower coordinates are valid
               stop("window area is zero for at least one observation window")
            }
            # Calculate the multivariate-normal normalisation integral for the current dimension
            normIntegral <- normIntegral * (pnorm((upperCoords[1:numObsWindows, dimIter] - sourceCoords[dimIter]) / normSD, 0.0, 1.0) - pnorm((lowerCoords[1:numObsWindows, dimIter] - sourceCoords[dimIter]) / normSD, 0.0, 1.0))
         }
      } else {
         # The observation windows are transects
         # Currently not supported
         stop("currently transect observation windows are not supported for multivariate normal source point models")
      }
      # Return the normalisation interal
      return(normIntegral)
   }
)

### 1.3. ==== Calculate whether a point within a set of observation windows ====
isWithinWindow <- nimbleFunction(
   run = function(
      curCoords = double(1),
      lowerCoords = double(2),
      upperCoords = double(2),
      tolDist = double(0),
      areAreas = double(0, default = 1)
   ) {
      ## 1.3.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(1))
      ## 1.3.2. Sanity test the inputs ----
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(lowerCoords)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      }
      # Assess the number of observation windows
      numObsWindows <- dim(lowerCoords)[1]
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      } else if(numObsWindows != dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates")
      }
      if(dimCoords != dim(upperCoords)[2]) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Ensure that the source coordinates have the correct dimensionality
      if(length(curCoords) != dimCoords) {
         stop("point coordinates do match the coordinates of the observation windows")
      }
      # Ensure that the tolerance distance is correctly set
      if(tolDist < 0.0) {
         stop("tolerance distance cannot be negative")
      }
      ## 1.3.3. Assess intersection between the points and the observation windows ----
      isInObsWindow <- numeric(length = numObsWindows, value = 1.0, recycle = TRUE)
      if(areAreas) {
         # The observation window is an area/volume so assess that it falls within
         # that volume
         for(dimIter in 1:dimCoords) {
            if(sum(lowerCoords[1:numObsWindows, dimIter] >= upperCoords[1:numObsWindows, dimIter]) > 0) {
               stop("upper coordinates must be greater than lower coordinates in each dimension")
            }
            isInObsWindow <- isInObsWindow * 
               numeric(value = rep(curCoords[dimIter], numObsWindows) >= lowerCoords[1:numObsWindows, dimIter] &
                   rep(curCoords[dimIter], numObsWindows) < upperCoords[1:numObsWindows, dimIter], length = numObsWindows)
         }
      } else {
         # The observation window is a line so test to see if the point falls on that
         # line
         isInObsWindow <- numeric(value = rep(curCoords[1], numObsWindows) >= lowerCoords[1:numObsWindows, 1] &
            rep(curCoords[1], numObsWindows) < upperCoords[1:numObsWindows, 1], length = numObsWindows)
         if(dimCoords > 1) {
            for(dimIter in 2:dimCoords) {
               isInObsWindow <- isInObsWindow * numeric(value = abs(rep(curCoords[dimIter], numObsWindows) - (
                  lowerCoords[1:numObsWindows, dimIter] + (upperCoords[1:numObsWindows, dimIter] - lowerCoords[1:numObsWindows, dimIter]) * (curCoords[1] - lowerCoords[1:numObsWindows, 1]) / (upperCoords[1:numObsWindows, 1] - lowerCoords[1:numObsWindows, 1])
               )) <= tolDist, length = numObsWindows)
            }
         }
      }
      return(isInObsWindow)
   }
)

### 1.4. ==== Find the nearest point in each observation window to a source point ====
nearestObsWindowPoint <- nimbleFunction(
   run = function(
      lowerCoords = double(2),
      upperCoords = double(2),
      sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
      areAreas = double(0, default = 1)
   ) {
      ## 1.4.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(2))
      ## 1.4.2. Sanity test the inputs ----
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(lowerCoords)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      }
      # Assess the number of observation windows
      numObsWindows <- dim(lowerCoords)[1]
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      } else if(numObsWindows != dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates")
      }
      if(dimCoords != dim(upperCoords)[2]) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Ensure that the source coordinates have the correct dimensionality
      if(length(sourceCoords) != dimCoords) {
         stop("source coordinates do match the coordinates of the observation windows")
      }
      ## 1.4.3. Calculate the set of nearest point to the source coordinates for each observation window ----
      nearCoords <- matrix(value = 0.0, nrow = numObsWindows, ncol = dimCoords, recycle = TRUE)
      if(areAreas) {
         # Observation windows are areas/volumes so calculate the nearest point appropriately
         for(dimIter in 1:dimCoords) {
            nearCoords[1:numObsWindows, dimIter] <- pmin(
               pmax(rep(sourceCoords[dimIter], numObsWindows), lowerCoords[1:numObsWindows, dimIter]),
               upperCoords[1:numObsWindows, dimIter])
         }
      } else {
         # Observation windows are transects so calculate the nearest point appropriately
         # Currently not supported
         stop("currently transect observation windows are not supported for multivariate normal source point models")
      }
      return(nearCoords)
   }
)

### 1.5. ==== Calculate the integral of the normal distribution function ====
# Utility function to calculate the integral of the normal distribution function
normalDistIntegral <- nimbleFunction(
   run = function(
      curVals = double(1)
   ) {
      ## 1.5.1. Specify the return type dimensionality ----
      returnType(double(1))
      ## 1.5.2. Calculate the integrated density ----
      outVals <- curVals * pnorm(curVals, 0.0, 1.0) + dnorm(curVals, 0.0, 1.0)
      return(outVals)
   }
)

### 1.6. ==== Calculate void probability: Binomial Point Pattern Parents Producing Poisson Pattern Children ====
# Function to calculate the probability of observing no points in a set of detection observation windows if
# the underlying process that generates the point pattern consists of the following two steps:
#     1. Parent points are generated according to an (inhomogenous) binomial point pattern process
#     2. Child points are produced according to an (inhomogenous) Poisson point pattern process that
#        is thinned according to an isotropic multivariate normal distribution.
# Only the child points are retained as part of the process.
voidProb_BinomialPPParent_PoisMNormSourcePPChild <- nimbleFunction(
   run = function(
      lowerHabitatWindow = double(2),
      upperHabitatWindow = double(2),
      habitatIntensity = double(1),
      lowerDetectWindow = double(2),
      upperDetectWindow = double(2),
      detectIntensity = double(1),
      detectSigma = double(0),
      areAreasHabitat = double(0, default = 1),
      areAreasDetect = double(0, default = 1),
      nHabWindowsN = double(0, default = -1),
      nDetectWindowsN = double(0, default = -1),
      numParents = double(0, default = 1)
   ) {
      ## 1.6.1. Specify the return type dimensionality ----
      returnType(double(0))
      ## 1.6.2. Sanity test the inputs ----
      # Check to ensure that the number of parents is correct
      if(numParents <= 0.0) {
         return(1.0)
      }
      # Retrieve the number of habitat windows
      numHabWindows <- trunc(nHabWindowsN)
      if(numHabWindows <= 0) {
         numHabWindows <- dim(lowerHabitatWindow)[1]
      }
      # Retrieve the number of dimensions in the coordinates
      dimCoords <- dim(lowerHabitatWindow)[2]
      if(numHabWindows <= 0) {
         # Use the dimensions of the lower habitat window matrix if the number of habitat windows
         # is set to below zero
         stop("invalid number of habitat-model observation windows")
      }
      if(dimCoords <= 0) {
         stop("invalid number of dimensions for the coordinates")
      }
      # Ensure that the dimensionality of the lower coordinates (habitat model) are correct
      if(dim(lowerHabitatWindow)[1] < numHabWindows) {
         stop("number of habitat-model observation windows not consistent between parameters")
      }
      # Ensure that the dimensionality of the upper coordinates (habitat model) are correct
      if(dim(upperHabitatWindow)[1] < numHabWindows) {
         stop("number of habitat-model observation windows not consistent between parameters")
      }
      if(dim(upperHabitatWindow)[2] != dimCoords) {
         stop("number of dimensions for the coordinates not consistent between parameters")
      }
      # Calculate the area of the habitat-model observation windows
      habitatWindowArea <- calcObsWindowSize(lowerHabitatWindow, upperHabitatWindow, areAreasHabitat)
      # Ensure that the intensity values for the habitat model are correct
      recHabitatIntensity <- numeric(value = habitatIntensity, length = numHabWindows, recycle = TRUE)
      recHabitatIntensity <- recHabitatIntensity * numeric(recHabitatIntensity > 0.0, length = numHabWindows, recycle = TRUE)
      if(sum(recHabitatIntensity) <= 0.0) {
         # Renormalise the habitat weights if all values are zero (this might not be desired behaviour)
         recHabitatIntensity <- rep(1.0, numHabWindows) / habitatWindowArea
      }
      # Ensure that the dimensionality of the detection windows (lower coordinates) is correct
      numDetectWindows <- trunc(nDetectWindowsN)
      if(numDetectWindows <= 0) {
         # Use the dimensions of the lower detection window matrix if the number of detection windows
         # is set to below zero
         numDetectWindows <- dim(lowerDetectWindow)[1]
      }
      if(numDetectWindows <= 0) {
         stop("invalid number of detection-model observation windows")
      }
      if(dim(lowerDetectWindow)[2] != dimCoords) {
         stop("number of dimensions for the coordinates not consistent between parameters")
      }
      # Ensure that the dimensionality of the detection windows (lower coordinates) is correct
      if(dim(lowerDetectWindow)[1] < numDetectWindows) {
         stop("number of detect-model observation windows not consistent between parameters")
      }
      # Ensure that the dimensionality of the detection windows (upper coordinates) is correct
      if(dim(upperDetectWindow)[1] < numDetectWindows) {
         stop("number of detect-model observation windows not consistent between parameters")
      }
      if(dim(upperDetectWindow)[2] != dimCoords) {
         stop("number of dimensions for the coordinates not consistent between parameters")
      }
      # Ensure that the intensity values for the detection model are correct
      recDetectIntensity <- numeric(value = detectIntensity, length = numDetectWindows, recycle = TRUE)
      recDetectIntensity <- recDetectIntensity * numeric(recDetectIntensity > 0.0, length = numDetectWindows, recycle = TRUE)
      if(sum(recDetectIntensity) <= 0.0) {
         # If the detection intensity is zero everywhere then the chance of not observing anything is one
         return(0.0)
      }
      # Ensure that the value for detection decay kernel
      if(detectSigma <= 0.0) {
         # These values should never occur but NIMBLE still sometimes proposes them so this instance needs to be handled
         # For the sake of these situations, we will treat zero (or lower) values of sigma as dirac delta filtering process
         # which mean that there will be no mass on points not directly on top of the source
         return(0.0)
      }
      ## 1.6.3.  Calculate the integral ----
      # Intialise an output value
      outVal <- -Inf
      # Initialise a temporary vector to hold integration over the habitat windows
      habWindowInt <- numeric(length = numHabWindows, value = 0.0)
      # Initialise a temporary vector to hold integration over the detection windows
      detectWindowInt <- numeric(length = numDetectWindows, value = 0.0)
      # Calculate the habitat normalisation value (integration of the parent point pattern intensity surface over the entire set of habitat windows)
      habNormVal <- sum(recHabitatIntensity * habitatWindowArea)
      if(areAreasHabitat & areAreasDetect) {
         for(detectIter in 1:numDetectWindows) {
            # For each detection window, integrate over the entire set of potential parent locations
            for(habIter in 1:numHabWindows) {
               habWindowInt[habIter] <- recHabitatIntensity[habIter] * prod(
                  normalDistIntegral((upperDetectWindow[detectIter, 1:dimCoords] - lowerHabitatWindow[habIter, 1:dimCoords]) / detectSigma) -
                  normalDistIntegral((upperDetectWindow[detectIter, 1:dimCoords] - upperHabitatWindow[habIter, 1:dimCoords]) / detectSigma) -
                  normalDistIntegral((lowerDetectWindow[detectIter, 1:dimCoords] - lowerHabitatWindow[habIter, 1:dimCoords]) / detectSigma) +
                  normalDistIntegral((lowerDetectWindow[detectIter, 1:dimCoords] - upperHabitatWindow[habIter, 1:dimCoords]) / detectSigma))
            }
            detectWindowInt[detectIter] <- recDetectIntensity[detectIter] * sum(habWindowInt[1:numHabWindows])
         }
         # Calculate the void probability for the composite intensity surface
         outVal <- pow(exp(-(
            # Calculate the full intensity integration for the composite intensity surface
            sum(detectWindowInt[1:numDetectWindows]) * pow(detectSigma, 2.0 * dimCoords) * pow(2.0 * pi, 0.5 * dimCoords) / habNormVal
         # Raise the calculated void probability to the power of each parent (to account for that no child of any parent was observed)
         )), trunc(numParents))
      } else {
         # Not currently able to do transect-based integration
         stop("transect-based integration not currently implemented")
      }
      return(outVal)
   }
)

### 1.7. ==== Lookup the index of an observation window that contains a set of source coordinates ====
# Function to input a matrix of source coordinates and return a vcetor of the indeces of the observation
# windows that the observation window falls within.
obsWindowIndex <- nimbleFunction(
   run = function(
      lowerCoords = double(2),
      upperCoords = double(2),
      sourceCoords = double(2)
   ) {
      ## 1.7.1. Specify the return type dimensionality ----
      returnType(double(1))
      ## 1.7.2. Sanity test the inputs ----
      # Ensure that the number of coordinates is valid
      dimCoords <- dim(lowerCoords)[2]
      if(dimCoords <= 0) {
         stop("invalid dimensionality of the coordinate system")
      }
      # Ensure that the number of windows is valid
      numWindows <- dim(lowerCoords)[1]
      if(numWindows <= 0) {
         stop("invalid number of observation windows")
      }
      # Ensure that the number of coordinates is consistent between the inputs
      if(dim(upperCoords)[2] != dimCoords | dim(sourceCoords)[2] != dimCoords) {
         stop("dimensionality of the inputs is not consisitent")
      }
      if(dim(upperCoords)[1] != numWindows) {
         stop("number of observation windows in the inputs is not consistent")
      }
      # Ensure that the number of points to get the indeces for is valid
      numPoints <- dim(sourceCoords)[1]
      if(numPoints <= 0) {
         stop("invalid number of source points")
      }
      ## 1.7.3. Calculate the indeces ----
      # Initialise a vector of output indeces
      outIndeces <- rep(0, numPoints)
      # Iterate over each of the source points
      for(pointIter in 1:numPoints) {
         # Initialise a flag denoting whether the corrct index has been found
         notFound <- 1
         # Iterate over the observation windows
         while(outIndeces[pointIter] < numWindows & notFound == 1) {
            outIndeces[pointIter] <- outIndeces[pointIter] + 1
            numAboveLowerCoords <- sum(sourceCoords[pointIter, 1:dimCoords] >= lowerCoords[outIndeces[pointIter], 1:dimCoords])
            numBelowUpperCoords <- sum(sourceCoords[pointIter, 1:dimCoords] < upperCoords[outIndeces[pointIter], 1:dimCoords])
            # Test to see if the source point falls within the current observation window
            if(numAboveLowerCoords == dimCoords & numBelowUpperCoords == dimCoords) {
               notFound <- 0
            }
         }
         if(notFound == 1) {
            # Index not found: the source point does not fall within an observation window
            outIndeces[pointIter] <- 0
         }
      }
      return(outIndeces)
   }
)

### 1.8. ==== Function to reduce the number of observation windows further than a given distance from a source coordinate ====
# Reduce the number of observation windows according to a distance criterion (can be used to speed up integration calculations)
obsWindowReduce <- nimbleFunction(
   run = function(
      lowerCoords = double(2),
      upperCoords = double(2),
      sourceCoords = double(1),
      maxDist = double(0),
      areAreas = double(0, default = 1)
   ) {
      ## 1.8.1. Specify the return type dimensionality ----
      returnType(double(2))
      ## 1.8.2. Sanity test the inputs ----
      # Ensure the number of coordinates
      dimCoords <- dim(lowerCoords)[2]
      if(dimCoords <= 0) {
         stop("invalid dimensionality of the coordinate system")
      }
      # Ensure that the number of windows is valid
      numWindows <- dim(lowerCoords)[1]
      if(numWindows <= 0) {
         stop("invalid number of observation windows")
      }
      # Ensure that the number of coordinates is consistent between the inputs
      if(dim(upperCoords)[2] != dimCoords | length(sourceCoords) != dimCoords) {
         stop("dimensionality of the inputs is not consisitent")
      }
      if(dim(upperCoords)[1] != numWindows) {
         stop("number of observation windows in the inputs is not consistent")
      }
      ## 1.8.3. Perform distance check ----
      # Initialise a usage vector
      useRow <- integer(value = c(1), length = numWindows)
      if(maxDist > 0.0) {
         # Calculate the distance from the source coordinate to the nearest point on the observation window
         nearestPoint <- nearestObsWindowPoint(lowerCoords, upperCoords, sourceCoords, areAreas)
         # Iterate over the nearest points and see if they fall within the Euclidean distance
         for(windowIter in 1:numWindows) {
            # Calculate the Euclidean distance
            distVal <- pow(sum(pow(sourceCoords[1:dimCoords] - nearestPoint[windowIter, 1:dimCoords], 2.0)), 0.5)
            # Check to see if the distance is less than the maximum-permissable distance
            if(distVal > maxDist) {
               useRow[windowIter] <- 0
            }
         }
      }
      ## 1.8.4. Create the reduced observation window matrix ----
      # Get the number of rows within the maximum specified distance
      numNewWindows <- sum(useRow)
      # Initalise an output matrix
      outputMat <- matrix(value = 0, nrow = numNewWindows, ncol = dimCoords * 2)
      newWindowIter <- 0
      # Fill the initialised matrix with the rows that fall within the designated distance
      for(windowIter in 1:numWindows) {
         if(useRow[windowIter]) {
            # If the row is to be used then copy the lower and upper coordinates into the relevant places
            # of the matrix
            newWindowIter <- newWindowIter + 1
            outputMat[newWindowIter, 1:dimCoords] <- lowerCoords[windowIter, 1:dimCoords]
            outputMat[newWindowIter, dimCoords + 1:dimCoords] <- upperCoords[windowIter, 1:dimCoords]
         }
      }
      # Return the output matrix
      return(outputMat)
   }
)

## 2. ------ DEFINE A NIMBLE CUSTOM DISTRIBUTION FOR THE (IN)HOMOGENOUS BINOMIAL POINT PROCESS ------

### 2.1. ==== Define the density function ====
# Define a function for the density function for the (in)homogenous binomial process
dbinomPP <- nimbleFunction(
   run = function(
      x = double(2),                               # Coordinate values to calculate the density
      numPoints = double(0),                       # The number of points in the binomial process
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      intensityWeights = double(1, default = 1),   # Values used in the intensity surface (by default a homogeous process is produced)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      log = integer(0, default = 0)                # If not 0 then return the log density
   ) {
      # Maximum allowable tolerance for testing to see if points fall on lines (only used when areAreas = 0)
      tolDist <- 6.661338e-16
      ## 2.1.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(0))
      ## 2.1.2. Sanity test the inputs ----
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(x)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      }
      # Ensure that the number of points corresponds to the dimensionality of the input coordinates
      if(dim(x)[1] < numPoints) {
         # If the number of rows in the output is less than the number of points simulated by the
         # binomial process: return a zero likelihood.  This is often not required as in nearly every
         # application of the binomial process these values will be the same but we need to handle the
         # occasional exception
         if(log) {
            return(-Inf)
         } else {
            return(0.0)
         }
      }
      if(numPoints == 0) {
         # If no points are simulated then the output is certain
         if(log) {
            return(0.0)
         } else {
            return(1.0)
         }
      } else if(numPoints < 0) {
         # Invalid number of points: set likelihood to zero
         if(log) {
            return(-Inf)
         } else {
            return(0.0)
         }
      }
      # Assess the number of observation windows
      numObsWindows <- trunc(numWindows)
      if(numObsWindows <= 0) {
         numObsWindows <- dim(lowerCoords)[1]
      }
      # Ensure that the number of observation windows is valid
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      }
      # Check that the number of observation windows is consistent (lower coordinates)
      if(numObsWindows > dim(lowerCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      # Check that the number of observation windows is consistent (upper coordinates)
      if(numObsWindows > dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      if(dim(lowerCoords)[2] != dimCoords | dim(upperCoords)[2] != dimCoords) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Calculate the area/length of the observation windows
      obsWindowSize <- calcObsWindowSize(
         matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         areAreas)
      # Recycle the intensity weights to match the number of observation windows
      recIntensityWeights <- numeric(length = numObsWindows, value = intensityWeights, recycle = TRUE)
      # Ensure that the intensity weights are valid
      if(sum(recIntensityWeights < 0.0) > 0) {
         # Invalid values for the intensity weights: set likelihood to zero
         if(log) {
            return(-Inf)
         } else {
            return(0.0)
         }
      }
      # Find the sum of the product of the intensity weights and area
      sumIntensity <- sum(recIntensityWeights * obsWindowSize)
      if(sumIntensity <= 0.0) {
         # Invalid values for the intensity weights: set likelihood to zero
         if(log) {
            return(-Inf)
         } else {
            return(0.0)
         }
      }
      ## 2.1.3. Calculate the likelihood of each point ----
      logPointDens <- rep(0.0, numPoints)
      for(pointIter in 1:numPoints) {
         curCoords <- x[pointIter, 1:dimCoords]
         # Retrieve the observation windows which the point falls within
         isInObsWindow <- isWithinWindow(curCoords,
            matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
            matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
            tolDist, areAreas)
         # Calculate the sum of the intensity
         pointSumIntensity <- sum(isInObsWindow * recIntensityWeights)
         if(pointSumIntensity <= 0.0) {
            if(log) {
               # Return the log scale density if requested
               return(-Inf)
            } else {
               # Return the natural scale density if requested
               return(0.0)
            }
         }
         logPointDens[pointIter] <- log(pointSumIntensity)
      }
      ## 2.1.4. Return the retrieved density ----
      outProb <- sum(logPointDens) - numPoints * log(sumIntensity)
      if(log == 0) {
         # Export the output probability on the natural scale if requested
         outProb <- exp(outProb)
      }
      return(outProb)
   }
)

### 2.2. ==== Define the sampling function ====
# Define a function to draw random coordinates from an (in)homogenous binomial process
rbinomPP <- nimbleFunction(
   run = function(
      n = integer(0),                              # Number of samples to draw from the distribution
      numPoints = double(0),                       # The number of points in the binomial process
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      intensityWeights = double(1, default = 1),   # Values used in the intensity surface (by default a homogeous process is produced)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numWindows = double(0, default = -1)         # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
   ) {
      ## 2.2.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(2))
      ## 2.2.2. Sanity test the inputs ----
      # Ensure that only one sample is requested
      if(n <= 0) {
         stop("the number of requested samples must be above zero")
      } else if(n > 1) {
         print("rbinomPP only allows n = 1; using n = 1")
      }
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(lowerCoords)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      }
      # Test the values for the number of points
      if(numPoints < 0) {
         stop("invalid number of points to simulate")
      } else if(numPoints == 0) {
         # Return an empty matrix (0 rows and dimCoords columns)
         return(matrix(nrow = 0, ncol = dimCoords, type = "double"))
      }
      # Asses the number of observation windows
      numObsWindows <- trunc(numWindows)
      if(numObsWindows <= 0) {
         numObsWindows <- dim(lowerCoords)[1]
      }
      # Ensure that the number of observation windows is valid
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      }
      # Check that the number of observation windows is consistent (lower coordinates)
      if(numObsWindows > dim(lowerCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      # Check that the number of observation windows is consistent (upper coordinates)
      if(numObsWindows > dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      if(dim(upperCoords)[2] != dimCoords) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Calculate the area/length of the observation windows
      obsWindowSize <- calcObsWindowSize(
         matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         areAreas)
      # Recycle the intensity weights to match the number of observation windows
      recIntensityWeights <- numeric(length = numObsWindows, value = intensityWeights, recycle = TRUE)
      # Ensure that the intensity weights are valid
      if(sum(recIntensityWeights < 0.0) > 0) {
         recIntensityWeights <- recIntensityWeights * numeric(value = recIntensityWeights > 0.0, length = numObsWindows)
      }
      # Weight the areas by the intensity weights
      areaIntensityWeights <- recIntensityWeights * obsWindowSize
      # Find the sum of the product of the intensity weights and area
      sumIntensity <- sum(areaIntensityWeights)
      if(sumIntensity <= 0.0) {
         # If the area weights sum to zero then instead set them all to their relative sizee
         areaIntensityWeights <- obsWindowSize
         sumIntensity <- sum(areaIntensityWeights)
      }
      ## 2.2.3. Generate random points ----
      # Generate observation window indeces for the output points
      # Currently this can't be done as a single call to rcat due to the way NIMBLE implements the categorical distribution
      obsWindowInd <- rep(0, numPoints)
      for(indIter in 1:numPoints) {
         obsWindowInd[indIter] <- rcat(1, areaIntensityWeights)
      }
      # Initialise an output matrix for the coordinates
      outCoordinates <- matrix(nrow = numPoints, ncol = dimCoords)
      if(areAreas){
         # The observation windows are areas/volumes so generate a set of random
         # coordinates within these volumes
         for(pointIter in 1:numPoints) {
            for(dimIter in 1:dimCoords) {
               outCoordinates[pointIter, dimIter] <- runif(1, min = lowerCoords[obsWindowInd[pointIter], dimIter], max = upperCoords[obsWindowInd[pointIter], dimIter])
            }
         }
      } else {
         # The observation windows are transects so generate a set of random coordinates
         # along their lengths
         for(pointIter in 1:numPoints) {
            coordProp <- runif(1, min = 0, max = 1)
            for(dimIter in 1:dimCoords) {
               outCoordinates[pointIter, dimIter] <- lowerCoords[obsWindowInd[pointIter], dimIter] + 
                  (upperCoords[obsWindowInd[pointIter], dimIter] - lowerCoords[obsWindowInd[pointIter], dimIter]) * coordProp
            }
         }
      }
      ## 2.2.4. Return the generated coordinates ----
      return(outCoordinates)
   }
)

### 2.3. ==== Define the density function for the single data-point version ====
dbinomPPSingle <- nimbleFunction(
	run = function(
		x = double(1),                               # Coordinate values to calculate the density
		lowerCoords = double(2),                     # The lower coordinate values of the observation windows
		upperCoords = double(2),                     # The upper coordinate values of the observation windows
		intensityWeights = double(1, default = 1),   # Values used in the intensity surface (by default a homogeous process is produced)
		areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
		numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
		log = integer(0, default = 0)                # If not 0 then return the log density
	) {
		## 2.3.1. Specify the return type dimensionality ----
		returnType(double(0))
		## 2.3.2. Create a temporary input matrix ----
		temporaryInput <- matrix(x, ncol = length(x), nrow = 1)
		## 2.3.3. Call the matrix-version of dbinomPP ----
		return(dbinomPP(temporaryInput, 1, lowerCoords, upperCoords, intensityWeights, areAreas, numWindows, log))
	}
)

### 2.4. ==== Define the sampling function for the single data-point version ====
rbinomPPSingle <- nimbleFunction(
	run = function(
		n = integer(0),                              # Number of samples to draw from the distribution
		lowerCoords = double(2),                     # The lower coordinate values of the observation windows
		upperCoords = double(2),                     # The upper coordinate values of the observation windows
		intensityWeights = double(1, default = 1),   # Values used in the intensity surface (by default a homogeous process is produced)
		areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
		numWindows = double(0, default = -1)         # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
	) {
		## 2.4.1. Specify the return dimensionality ----
		returnType(double(1))
		## 2.4.2. Create a temporary output matrix ----
		# Retrieve the number of dimensions
		dimCoords <- dim(lowerCoords)[2]
		# Call the matrix-version of rbinomPP
		temporaryOutput <- rbinomPP(n, 1, lowerCoords, upperCoords, intensityWeights, areAreas, numWindows)
		## 2.4.3. Return a slice of the output matrix ----
		return(temporaryOutput[1, 1:dimCoords])
	}
)

## 3. ------ DEFINE A NIMBLE CUSTOM DISTRIBUTION FOR THE MULTIVARIATE-NORMAL SOURCE POINT INHOMOGENOUS BINOMIAL POINT PROCESS (AND DERIVATIVES) ------

### 3.1. ==== Define the density function ====
# Define a function for the density function for the multivariate-normal source point inhomogenous binomial process
dbinomMNormSourcePP <- nimbleFunction(
   run = function(
      x = double(2),                               # Coordinate values to calculate the density
      numPoints = double(0),                       # The number of points in the binomial process
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      localEvalParam = double(0, default = -1),    # Parameter that controls the maximum distance that an observation window can be from x to be considered as a possible destination (-1 is interpredted as +Inf)
      log = integer(0, default = 0)                # If not 0 then return the log density
   ) {
      # Maximum allowable tolerance for testing to see if points fall on lines (only used when areAreas = 0)
      tolDist <- 6.661338e-16
      ## 3.1.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(0))
      ## 3.1.2. Sanity test the inputs ----
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(x)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      } else if(dimCoords != 2) {
         stop("currently the bivariate-normal source point inhomogenous point process model is only defined for 2-dimensional processes")
      }
      # Ensure that the number of points corresponds to the dimensionality of the input coordinates
      if(dim(x)[1] < numPoints) {
         # If the number of rows in the output is lower than the number of points simulated by the
         # binomial process: return a zero likelihood.  This is often not required as in nearly every
         # application of the binomial process these values will be the same but we need to handle the
         # occasional
         if(log) {
            # Return the log scale density if requested
            return(-Inf)
         } else {
            # Return the natural scale density if requested
            return(0.0)
         }
      }
      if(numPoints == 0) {
         # If no points are simulated then the output is certain
         if(log) {
            # Return the log scale density if requested
            return(0.0)
         } else {
            # Return the natural scale density if requested
            return(1.0)
         }
      } else if(numPoints < 0) {
         # Invalid number of points: set likelihood to zero
         if(log) {
            # Return the log scale density if requested
            return(-Inf)
         } else {
            # Return the natural scale density if requested
            return(0.0)
         }
      }
      # Asses the number of observation windows
      numObsWindows <- trunc(numWindows)
      if(numObsWindows <= 0) {
         numObsWindows <- dim(lowerCoords)[1]
      }
      # Ensure that the number of observation windows is valid
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      }
      # Check that the number of observation windows is consistent (lower coordinates)
      if(numObsWindows > dim(lowerCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      # Check that the number of observation windows is consistent (upper coordinates)
      if(numObsWindows > dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      # Ensure that the lower and upper coordinates have the correct dimension structure
      if(dim(lowerCoords)[2] != dimCoords | dim(upperCoords)[2] != dimCoords) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Perform the window reduction according to the local evaluation criterion
      redMatrix <- obsWindowReduce(lowerCoords, upperCoords, sourceCoords, localEvalParam, areAreas)
      numObsWindows <- dim(redMatrix)[1]
      if(numObsWindows <= 0) {
         # If there are no observation windows after the reduction process then return a zero likelihood
         if(log) {
            # Return the log scale density if requested
            return(-Inf)
         } else {
            # Return the natural scale density if requested
            return(0.0)
         }
      }
      inLowerCoords <- matrix(nrow = numObsWindows, ncol = dimCoords)
      inUpperCoords <- matrix(nrow = numObsWindows, ncol = dimCoords)
      inLowerCoords[1:numObsWindows, 1:dimCoords] <- redMatrix[1:numObsWindows, 1:dimCoords]                              # Retrieve the lower coordinates of the reduced windows
      inUpperCoords[1:numObsWindows, 1:dimCoords] <- redMatrix[1:numObsWindows, (dimCoords + 1):(dimCoords + dimCoords)]  # Retrieve the upper coordinates of the reduced windows
      # Recycle the intensity weights to match the number of observation windows
      recIntensityWeights <- numeric(length = numObsWindows, value = intensityWeights, recycle = TRUE)
      # Ensure that the intensity weights are valid
      if(sum(recIntensityWeights < 0.0) > 0) {
         # Invalid values for the intensity weights: set likelihood to zero
         if(log) {
            # Return the log scale density if requested
            return(-Inf)
         } else {
            # Return the natural scale density if requested
            return(0.0)
         }
      }
      # Test the validity of the source coordinates
      if(length(sourceCoords) != dimCoords) {
         stop("dimensioality of source coordinates are not consistent")
      }
      # Test the validiity of the decay parameter
      if(normSD <= 0.0) {
         if(log) {
            # Return the log scale density if requested
            return(-Inf)
         } else {
            # Return the natural scale density if requested
            return(0.0)
         }
      }
      ## 3.1.3. Calculate the likelihood of each point ----
      # Calculate the integration of the decay kernel (weighted by intensity) for each observation window
      obsWindowNorm <- calcMNormSourceInt(
         matrix(inLowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         matrix(inUpperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         sourceCoords, normSD, areAreas) * recIntensityWeights
      totalNorm <- sum(obsWindowNorm)
      logPointDens <- rep(0.0, numPoints)
      if(totalNorm <= 0.0) {
         # Return a zero likelihood if the intensity values sum to zero
         if(log) {
            # Return the log scale density if requested
            return(-Inf)
         } else {
            # Return the natural scale density if requested
            return(0.0)
         }
      } else {
         # Calculate the likelihood for each point
         for(pointIter in 1:numPoints) {
            curCoords <- x[pointIter, 1:dimCoords]
            # Retrieve the observation windows which the point falls within
            isInObsWindow <- isWithinWindow(curCoords,
               matrix(inLowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
               matrix(inUpperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
               tolDist, areAreas)
            # Calculate the distance decay component of the current point from the source (the exponent
            # of the isotropic multivariate normal diistribution)
            distDecay <- 0.0
            for(dimIter in 1:dimCoords) {
               distDecay <- distDecay + (curCoords[dimIter] - sourceCoords[dimIter]) * (curCoords[dimIter] - sourceCoords[dimIter])
            }
            distDecay <- distDecay * (-0.5 / (normSD * normSD))
            # Calculate the sum of the intensity for each point
            pointSumIntensity <- sum(isInObsWindow * recIntensityWeights * exp(distDecay))
            if(pointSumIntensity <= 0.0) {
               if(log) {
                  # Return the log scale density if requested
                  return(-Inf)
               } else {
                  # Return the natural scale density if requested
                  return(0.0)
               }
            }
            logPointDens[pointIter] <- log(pointSumIntensity)
         }
      }
      ## 3.1.4. Return the retrieved density ----
      outProb <- sum(logPointDens) - numPoints * log(totalNorm)
      if(log == 0) {
         # Export the output probability on the natural scale if requested
         outProb <- exp(outProb)
      }
      return(outProb)
   }
)

### 3.2. ==== Define the sampling function ====
# Define a function to draw random coordinate from the multivariate-normal source point inhomogenous binomial process
rbinomMNormSourcePP <- nimbleFunction(
   run = function(
      n = integer(0),                              # Number of samples to draw from the distribution
      numPoints = double(0),                       # The number of points in the binomial process
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      localEvalParam = double(0, default = -1)     # Parameter that controls the maximum distance that an observation window can be from x to be considered as a possible destination (-1 is interpredted as +Inf)
   ) {
      ## 3.2.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(2))
      ## 3.2.2. Sanity test the inputs ----
      # Ensure that only one sample is requested
      if(n <= 0) {
         stop("the number of requested samples must be above zero")
      } else if(n > 1) {
         print("rbinomMNormSourcePP only allows n = 1; using n = 1")
      }
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(lowerCoords)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      }
      # Test the values for the number of points
      if(numPoints < 0) {
         stop("invalid number of points to simulate")
      } else if(numPoints == 0) {
         # Return an empty matrix (0 rows and dimCoords columns)
         return(matrix(nrow = 0, ncol = dimCoords, type = "double"))
      }
      # Assess the number of observation windows
      numObsWindows <- trunc(numWindows)
      if(numObsWindows <= 0) {
         numObsWindows <- dim(lowerCoords)[1]
      }
      # Ensure that the number of observation windows is valid
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      }
      # Check that the number of observation windows is consistent (lower coordinates)
      if(numObsWindows > dim(lowerCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      # Check that the number of observation windows is consistent (upper coordinates)
      if(numObsWindows > dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      if(dim(upperCoords)[2] != dimCoords) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Recycle the intensity weights to match the number of observation windows
      recIntensityWeights <- numeric(length = numObsWindows, value = intensityWeights, recycle = TRUE)
      # Ensure that the intensity weights are valid
      if(sum(recIntensityWeights < 0.0) > 0) {
         recIntensityWeights <- recIntensityWeights * numeric(value = recIntensityWeights > 0.0, length = numObsWindows)
      }
      # Test the validity of the source coordinates
      if(length(sourceCoords) != dimCoords) {
         stop("dimensionality of source coordinates are not consistent")
      }
      # Test the validiity of the decay parameter
      if(normSD <= 0.0) {
         stop("invalid values for the decay parameter")
      }
      # Print error message if the local evaluation parameters are used in the simulation function
      if(localEvalParam > 0.0) {
         print("local evaluation parameters not used in the simulation function")
      }
      ## 3.2.3. Generate random points ----
      # Current a stratified rejection sampling approach is employed to generate random variables from this process.
      # This can become inefficient if the observation windows because very large compared to the range of the decay
      # kernel.  This could be optimised better in future versions.
      obsWindowInd <- rep(1, numPoints)
      if(numObsWindows > 1) {
         # Weight the observation windows by their intensity weights multiplied by the integral of the decay function
         integratedSource <- calcMNormSourceInt(
            matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
            matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
            sourceCoords, normSD, areAreas)
         areaIntensityWeights <- recIntensityWeights * integratedSource
         # Find the sum of the product of the intensity weights and area
         sumIntensity <- sum(areaIntensityWeights)
         if(sumIntensity <= 0.0) {
            areaIntensityWeights <- integratedSource
            sumIntensity <- sum(areaIntensityWeights)
         }
         # Generate observation window indeces for the output
         # Currently this can't be done as a single call to rcat due to the way NIMBLE implements the categorical distribution
         for(indIter in 1:numPoints) {
            obsWindowInd[indIter] <- rcat(1, areaIntensityWeights)
         }
      }
      # Initialise an output matrix for the coordinates
      outCoordinates <- matrix(nrow = numPoints, ncol = dimCoords)
      # Calculate the nearest point to the source coordinates in each observation window
      nearestPoints <- nearestObsWindowPoint(
         matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         sourceCoords, areAreas)
      # Calculate the value of the decay kernel at the nearest point in each observation window
      decayNearest <- rep(0.0, numObsWindows)
      for(dimIter in 1:dimCoords) {
         decayNearest <- decayNearest + (nearestPoints[1:numObsWindows, dimIter] - sourceCoords[dimIter]) * (nearestPoints[1:numObsWindows, dimIter] - sourceCoords[dimIter])
      }
      decayNearest <- (decayNearest * -0.5) / (normSD * normSD)
      decayNearest <- exp(decayNearest)
      if(areAreas) {
         # The observation windows are areas/volumes so generate a set of random
         # coordinates within these volumes
         for(pointIter in 1:numPoints) {
            # Create a set of test coordinates
            testCoords <- runif(dimCoords, min = lowerCoords[obsWindowInd[pointIter], 1:dimCoords], max = upperCoords[obsWindowInd[pointIter], 1:dimCoords])
            # Assess the value of the decay parameter
            decayTest <- (testCoords - sourceCoords) * (testCoords - sourceCoords)
            decayVal <- exp((-0.5 * sum(decayTest)) / (normSD * normSD))
            randVal <- runif(1, min = 0, max = 1)
            while(randVal > (decayVal / decayNearest[obsWindowInd[pointIter]])) {
               # Create a set of test coordinates
               testCoords <- runif(dimCoords, min = lowerCoords[obsWindowInd[pointIter], 1:dimCoords], max = upperCoords[obsWindowInd[pointIter], 1:dimCoords])
               # Assess the value of the decay parameter
               decayTest <- (testCoords - sourceCoords) * (testCoords - sourceCoords)
               decayVal <- exp((-0.5 * sum(decayTest)) / (normSD * normSD))
               randVal <- runif(1, min = 0, max = 1)
            }
            # Save the test coordinates once they have been accepted
            outCoordinates[pointIter, 1:dimCoords] <- testCoords
         }
      } else {
         # The observation windows are transects so generate a set of random coordinates
         # along their lengths
         for(pointIter in 1:numPoints) {
            # Create a set of test coordinates
            testProp <- runif(1, min = 0, max = 1)
            testCoords <- lowerCoords[obsWindowInd[pointIter], 1:dimCoords] + (upperCoords[obsWindowInd[pointIter], 1:dimCoords] - lowerCoords[obsWindowInd[pointIter], 1:dimCoords]) * testProp
            # Assess the value of the decay parameter
            decayTest <- (testCoords - sourceCoords) * (testCoords - sourceCoords)
            decayVal <- exp((-0.5 * sum(decayTest)) / (normSD * normSD))
            randVal <- runif(1, min = 0, max = 1)
            while(randVal > (decayVal / decayNearest[obsWindowInd[pointIter]])) {
               # Create a set of test coordinates
               testProp <- runif(1, min = 0, max = 1)
               testCoords <- lowerCoords[obsWindowInd[pointIter], 1:dimCoords] + (upperCoords[obsWindowInd[pointIter], 1:dimCoords] - lowerCoords[obsWindowInd[pointIter], 1:dimCoords]) * testProp
               # Assess the value of the decay parameter
               decayTest <- (testCoords - sourceCoords) * (testCoords - sourceCoords)
               decayVal <- exp((-0.5 * sum(decayTest)) / (normSD * normSD))
               randVal <- runif(1, min = 0, max = 1)
            }
            # Save the test coordinates once they have been accepted
            outCoordinates[pointIter, 1:dimCoords] <- testCoords
         }
      }
      ## 3.2.4. Return the generated coordinates ----
      return(outCoordinates)
   }
)

### 3.3. ==== Define the density function for the single data-point version ====
dbinomMNormSourcePPSingle <- nimbleFunction(
   run = function(
      x = double(1),                               # Coordinate values to calculate the density
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      localEvalParam = double(0, default = -1),    # Parameter that controls the maximum distance that an observation window can be from x to be considered as a possible destination (-1 is interpredted as +Inf)
      log = integer(0, default = 0)                # If not 0 then return the log density
   ) {
		## 3.3.1. Specify the return type dimensionality ----
      returnType(double(0))
      ## 3.3.2. Create a temporary input matrix ----
      temporaryInput <- matrix(x, ncol = length(x), nrow = 1)
      ## 3.3.3. Call the matrix-version of dbinomMNormSourcePP ----
      return(dbinomMNormSourcePP(temporaryInput, 1, lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, localEvalParam, log))
   }
)

### 3.4. ==== Define the sampling function for the single data-point version ====
rbinomMNormSourcePPSingle <- nimbleFunction(
   run = function(
      n = integer(0),                              # Number of samples to draw from the distribution
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      localEvalParam = double(0, default = -1)     # Parameter that controls the maximum distance that an observation window can be from x to be considered as a possible destination (-1 is interpredted as +Inf)
   ) {
      ## 3.4.1. Specify the return type dimensionality ----
      returnType(double(1))
      ## 3.4.2. Sanity test the inputs ----
      # Ensure that only one sample is requested
      if(n <= 0) {
         stop("the number of requested samples must be above zero")
      } else if(n > 1) {
         print("rbinomMNormSourcePPSingle only allows n = 1; using n = 1")
      }
      ## 3.4.3. Create a temporary output matrix ----
      # Retrieve the number of coordinates
      dimCoords <- dim(lowerCoords)[2]
      # Call the matrix-version of rbinomMNormSourcePP
      temporaryOutput <- rbinomMNormSourcePP(1, 1, lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, localEvalParam)
      ## 3.4.3. Return a slice of the output matrix ----
      return(temporaryOutput[1, 1:dimCoords])
   }
)

### 3.5. ==== Define the density function for the vectorised source point version ====
dbinomMNormSourcePPMulti <- nimbleFunction(
   run = function(
      x = double(2),                               # Coordinate values to calculate the density
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      sourceCoords = double(2),                    # The coordinates of the source location (the origin of the decay kernel) for each individual
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      localEvalParam = double(0, default = -1),    # Parameter that controls the maximum distance that an observation window can be from x to be considered as a possible destination (-1 is interpredted as +Inf)
      log = integer(0, default = 0)                # If not 0 then return the log density
   ) {
      ## 3.5.1. Specify the return type dimensionality ----
      returnType(double(0))
      ## 3.5.2. Sanity test the inputs ----
      # Retrieve the number of dimensions of the point coordinates
      dimCoords <- dim(x)[2]
      # Retrieve the number of points to calculate the AC displacement probability for
      numPoints <- dim(x)[1]
      if(dim(sourceCoords)[1] != numPoints) {
         stop("inconsistent number of rows between the distribution target and the source coordinates")
      }
      if(dim(sourceCoords)[2] != dimCoords) {
         stop("inconsistent number of columns  between the distribution target and the source coordinates")
      }
      # Ensure that the number of points is valid
      if(numPoints <= 0) {
         stop("invalid number of points")
      }
      # Ensure that the number of dimensions is valid
      if(dimCoords <= 0) {
         stop("invalid number of dimensions for the point coordinate system")
      }
      ## 3.5.3. Calculate the likelihood for each displacement ----
      # Initialise a vector to hold the log-probabilities for each displacement
      outProbVec <- numeric(length = numPoints, value = 0.0)
      for(pointIter in 1:numPoints) {
         # Calculate the log-likelihood of the displacement
         outProbVec[pointIter] <- dbinomMNormSourcePPSingle(x[pointIter, 1:dimCoords], lowerCoords, upperCoords, sourceCoords[pointIter, 1:dimCoords], normSD, intensityWeights, areAreas, numWindows, localEvalParam, 1)
      }
      ## 3.5.4. Aggregate into a single likelihood ----
      outProb <- sum(outProbVec[1:numPoints])
      if(log == 0) {
         # Output the likelihood on the natural scale if requested
         outProb <- exp(outProb)
      }
      return(outProb)
   }
)

### 3.6. ==== Define the source function for the vectorised source point version ====
rbinomMNormSourcePPMulti <- nimbleFunction(
   run = function(
      n = integer(0),                              # Number of samples to draw from the distribution
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      sourceCoords = double(2),                    # The coordinates of the source location (the origin of the decay kernel) for each individual
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      localEvalParam = double(0, default = -1)     # Parameter that controls the maximum distance that an observation window can be from x to be considered as a possible destination (-1 is interpredted as +Inf)
   ) {
      ## 3.6.1. Specify the return type dimensionality ----
      returnType(double(2))
      ## 3.6.2. Sanity test the inputs ----
      # Ensure that only one sample is requested
      if(n <= 0) {
         stop("the number of requested samples must be above zero")
      } else if(n > 1) {
         print("rbinomMNormSourcePPMulti only allows n = 1; using n = 1")
      }
      # Retrieve the number of points to generate
      numPoints <- dim(sourceCoords)[1]
      # Retrieve the number of dimensions
      dimCoords <- dim(sourceCoords)[2]
      # Ensure that the number of points is valid
      if(numPoints <= 0) {
         stop("invalid number of points")
      }
      # Ensure that the number of dimensions is valid
      if(dimCoords <= 0) {
         stop("invalid number of dimensions for the point coordinate system")
      }
      ## 3.6.3. Generate the displacements ----
      # Initialise a matrix to hold the outputs
      outPoints <- matrix(value = 0, nrow = numPoints, ncol = dimCoords)
      for(pointIter in 1:numPoints) {
         # Generate output points from successive calls to the single-source version
         outPoints[pointIter, 1:dimCoords] <- rbinomMNormSourcePPSingle(1, lowerCoords, upperCoords, sourceCoords[pointIter, 1:dimCoords], normSD, intensityWeights, areAreas, numWindows, localEvalParam)
      }
      ## 3.6.4. Return the generated output coordinates ----
      return(outPoints)
   }
)

### 3.7. ==== Define a desnsity function for the binomial point process model (single point, mixture model) ====
dbinomMNormSourcePP_dbinomPP_SingleMixture <- nimbleFunction(
   run = function(
      x = double(1),                               # Coordinate values to calculate the density
      mixtureParam = double(0),                    # Parameter to control the mixture frequency (dbinomPP freqeuncy)
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      localEvalParam = double(0, default = -1),    # Parameter that controls the maximum distance that an observation window can be from x to be considered as a possible destination (-1 is interpredted as +Inf)
      log = integer(0, default = 0)                # If not 0 then return the log density
   ) {
      ## 3.7.1. Define the return type ----
      returnType(double(0))
      ## 3.7.2. Sanity check the inputs ----
      # Ensure that the mixture parameter is between 0 and 1
      if(mixtureParam < 0.0 | mixtureParam > 1.0) {
         # Return zero likelihood if it is not
         if(log) {
            return(-Inf)
         } else {
            return(0.0)
         }
      }
      ## 3.7.3. Calculate the likelihood based on the weighted densities ----
      # Weight the densities based on the mixture parameter
      outProb <- mixtureParam * dbinomPPSingle(x, lowerCoords, upperCoords, intensityWeights, areAreas, numWindows, 0) +
         (1.0 - mixtureParam) * dbinomMNormSourcePPSingle(x, lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, 0)
      if(log) {
         outProb <- log(outProb)
      }
      return(outProb)
   }
)

### 3.8. ==== Define a function to draw coordinates from a binomial point process model (single point, mixture model) ====
rbinomMNormSourcePP_dbinomPP_SingleMixture <- nimbleFunction(
   run = function(
      n = integer(0),                              # Number of samples to draw from the distribution
      mixtureParam = double(0),                    # Parameter to control the mixture frequency (dbinomPP freqeuncy)      
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      localEvalParam = double(0, default = -1)     # Parameter that controls the maximum distance that an observation window can be from x to be considered as a possible destination (-1 is interpredted as +Inf)
   ) {
      ## 3.7.1. Define the return type ----
      returnType(double(1))
      ## 3.7.2. Determine which distribution to draw from ----
      # Initialise an output set of coordinates
      outCoords <- numeric(length = length(sourceCoords))
      if(runif(1, 0, 1) < mixtureParam) {
         # Draw from the inhomogenous binomial point process
         outCoords <- rbinomPPSingle(n, lowerCoords, upperCoords, intensityWeights, areAreas, numWindows)
      } else {
         # Draw from the multivariate-normal thinned binomial point process
         outCoords <- rbinomMNormSourcePPSingle(n, lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, localEvalParam)
      }
      return(outCoords)
   }
)

## 4. ------ DEFINE A NIMBLE CUSTOM DISTRIBUTION FOR THE (IN)HOMOGENOUS POISSON POINT PROCESS ------

### 4.1. ==== Define the density function ====
# Define a function for the density function for the (in)homogenous Poisson point process
dpoisPP <- nimbleFunction(
   run = function(
      x = double(2),                               # Coordinate values to calculate the density
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      intensityWeights = double(1, default = 1),   # Values used in the intensity surface (by default a homogeous process is produced)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numSamples = double(0, default = -1),        # Number of points (if negative then the number of rows in x is used to define this value)
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      allowZero = double(0, default = 1),          # Allow zero samples? Non-zero values denote that there can be
      log = integer(0, default = 0)                # If not 0 then return the log density
   ) {
      # Maximum allowable tolerance for testing to see if points fall on lines (only used when areAreas = 0)
      tolDist <- 6.661338e-16
      ## 4.1.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(0))
      ## 4.1.2. Sanity test the inputs ----
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(x)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      }
      # Ensure that the number of points corresponds to the dimensionality of the input coordinates
      numPoints <- trunc(numSamples)
      if(numPoints < 0) {
         numPoints <- dim(x)[1]
      }
      if(numPoints < 0) {
         # This code should never be called it is a hangover from a template function
         # Invalid number of points: set likelihood to zero
         if(log) {
            return(-Inf)
         } else {
            return(0.0)
         }
      }
      # Assess the number of observation windows
      numObsWindows <- trunc(numWindows)
      if(numObsWindows <= 0) {
         numObsWindows <- dim(lowerCoords)[1]
      }
      # Ensure that the number of observation windows is valid
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      }
      # Check that the number of observation windows is consistent (lower coordinates)
      if(numObsWindows > dim(lowerCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      # Check that the number of observation windows is consistent (upper coordinates)
      if(numObsWindows > dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      if(dim(lowerCoords)[2] != dimCoords | dim(upperCoords)[2] != dimCoords) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Calculate the area/length of the observation windows
      obsWindowSize <- calcObsWindowSize(
         matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         areAreas)
      # Recycle the intensity weights to match the number of observation windows
      recIntensityWeights <- numeric(length = numObsWindows, value = intensityWeights, recycle = TRUE)
      # Ensure that the intensity weights are valid
      if(sum(recIntensityWeights < 0.0) > 0) {
         # Invalid values for the intensity weights: set likelihood to zero
         if(log) {
            return(-Inf)
         } else {
            return(0.0)
         }
      }
      # Find the sum of the product of the intensity weights and area
      sumIntensity <- sum(recIntensityWeights * obsWindowSize)
      if(sumIntensity <= 0.0) {
         # If there are no points then this outcome is certain when the sum of the intensity
         # surface is zero
         if(numPoints == 0) {
            if(log) {
               return(0.0)
            } else {
               return(1.0)
            }
         # If there is at least one points then the outcome is impossible when the sum of the
         # intensity surface is zero
         } else {
            if(log) {
               return(-Inf)
            } else {
               return(0.0)
            }
         }
      }
      ## 4.1.3. Calculate the likelihood of each point ----
      logPointDens <- rep(0.0, numPoints)
      outProb <- -sumIntensity
      if(numPoints > 0) {
         for(pointIter in 1:numPoints) {
            curCoords <- x[pointIter, 1:dimCoords]
            # Retrieve the observation windows which the point falls within
            isInObsWindow <- isWithinWindow(curCoords,
               matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
               matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
               tolDist, areAreas)
            # Calculate the sum of the intensity
            pointSumIntensity <- sum(isInObsWindow * recIntensityWeights)
            if(pointSumIntensity <= 0.0) {
               if(log) {
                  # Return the log scale density if requested
                  return(-Inf)
               } else {
                  # Return the natural scale density if requested
                  return(0.0)
               }
            }
            logPointDens[pointIter] <- log(pointSumIntensity)
         }
         outProb <- sum(logPointDens) - sumIntensity# - sum(log(1:numPoints))
      } else if(allowZero == 0.0) {
         # Zero values have no likelihood weight
         if(log) {
            return(-Inf)
         } else {
            return(0.0)
         }
      }
      ## 4.1.4. Correct for truncation (if appropriate) ----
      if(allowZero == 0.0) {
         # Account for zero truncation
         outProb <- outProb - log(1.0 - exp(-sumIntensity))
      }
      ## 4.1.5. Return the retrieved density ----
      if(log == 0) {
         # Export the output probability on the natural scale if requested
         outProb <- exp(outProb)
      }
      return(outProb)
   }
)

### 4.2. ==== Define the sampling function ====
# Define a function to draw random coordinates from an (in)homogenous Poisson point process
rpoisPP <- nimbleFunction(
   run = function(
      n = integer(0),                              # Number of samples to draw from the distribution
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      intensityWeights = double(1, default = 1),   # Values used in the intensity surface (by default a homogeous process is produced)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numSamples = double(0, default = -1),        # Number of points (if negative then the number of rows in x is used to define this value)
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      allowZero = double(0, default = 1)           # Allow zero samples? Non-zero values denote that there can be
   ) {
      ## 4.2.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(2))
      ## 4.2.2. Sanity test the inputs ----
      # Ensure that only one sample is requested
      if(n <= 0) {
         stop("the number of requested samples must be above zero")
      } else if(n > 1) {
         print("rpoisPP only allows n = 1; using n = 1")
      }
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(lowerCoords)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      }
      # Assess the number of observation windows
      numObsWindows <- trunc(numWindows)
      if(numObsWindows <= 0) {
         numObsWindows <- dim(lowerCoords)[1]
      }
      # Ensure that the number of observation windows is valid
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      }
      # Check that the number of observation windows is consistent (lower coordinates)
      if(numObsWindows > dim(lowerCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      # Check that the number of observation windows is consistent (upper coordinates)
      if(numObsWindows > dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      if(dim(upperCoords)[2] != dimCoords) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Calculate the area/length of the observation windows
      obsWindowSize <- calcObsWindowSize(
         matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         areAreas)
      # Recycle the intensity weights to match the number of observation windows
      recIntensityWeights <- numeric(length = numObsWindows, value = intensityWeights, recycle = TRUE)
      # Ensure that the intensity weights are valid
      if(sum(recIntensityWeights < 0.0) > 0) {
         recIntensityWeights <- recIntensityWeights * numeric(value = recIntensityWeights > 0.0, length = numObsWindows)
      }
      # Weight the areas by the intensity weights
      areaIntensityWeights <- recIntensityWeights * obsWindowSize
      sumAreaIntensityWeights <- sum(areaIntensityWeights)
      if(sumAreaIntensityWeights <= 0.0) {
         return(matrix(nrow = 0, ncol = dimCoords, type = "double"))
      }
      ## 4.2.3. Generate a random number of points ----
      # Over-ride the random number of points if the number of samples has been set
      numPoints <- trunc(numSamples)
      if(numSamples < 0.0) {
         # If the number of samples has not been set then draw it from a Poisson distribution
         numPoints <- rpois(1, sumAreaIntensityWeights)
         if(allowZero == 0.0) {
            # If zeros are not possible then keep redrawing samples until a non-zero sample is drawn
            while(numPoints == 0) {
               numPoints <- rpois(1, sumAreaIntensityWeights)
            }
         }
      }
      if(numPoints < 0) {
         stop("invalid number of points to simulate")
      } else if(numPoints == 0) {
         if(allowZero == 0.0) {
            # Incompatible parameterisation
            stop("zero number of points requested for truncated distribution that excludes zero")
         }
         # Return an empty matrix (0 rows and dimCoords columns)
         return(matrix(nrow = 0, ncol = dimCoords, type = "double"))
      }
      ## 4.2.4. Generate random locations for the points ----
      # Sample the points according to the binomial point-process model
      outCoordinates <- rbinomPP(1, numPoints,
         matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         recIntensityWeights, areAreas)
      return(outCoordinates)
   }
)

## 5. ------ DEFINE A NUMBLE CUSTOM DISTRIBUTION FOR THE MULTIVARIATE-NORMAL SOURCE POINT INHOMOGENOUS POISSON POINT PROCESS ------

### 5.1. ==== Define the density function ====
# Define a function for the density function for the multivariate-normal source point inhomogenous Poisson point process
dpoisMNormSourcePP <- nimbleFunction(
   run = function(
      x = double(2),                               # Coordinate values to calculate the density
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numSamples = double(0, default = -1),        # Number of points (if negative then the number of rows in x is used to define this value)
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      allowZero = double(0, default = 1),          # Allow zero samples? Non-zero values denote that there can be
      localEvalParam = double(0, default = -1),    # Parameter that controls the maximum distance that an observation window can be from x to be considered as a possible destination (-1 is interpredted as +Inf)
      log = integer(0, default = 0)                # If not 0 then return the log density
   ) {
      # Maximum allowable tolerance for testing to see if points fall on lines (only used when areAreas = 0)
      tolDist <- 6.661338e-16
      ## 5.1.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(0))
      ## 5.1.2. Sanity test the inputs ----
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(x)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      } else if(dimCoords != 2) {
         stop("currently the bivariate-normal source point inhomogenous point process model is only defined for 2-dimensional processes")
      }
      # Ensure that the number of points corresponds to the dimensionality of the input coordinates
      numPoints <- trunc(numSamples)
      if(numPoints < 0) {
         numPoints <- dim(x)[1]
      }
      if(numPoints < 0) {
         # Invalid number of points: set likelihood to zero
         if(log) {
            # Return the log scale density if requested
            return(-Inf)
         } else {
            # Return the natural scale density if requested
            return(0.0)
         }
      }
      # Assess the number of observation windows
      numObsWindows <- trunc(numWindows)
      if(numObsWindows <= 0) {
         numObsWindows <- dim(lowerCoords)[1]
      }
      # Ensure that the number of observation windows is valid
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      }
      # Check that the number of observation windows is consistent (lower coordinates)
      if(numObsWindows > dim(lowerCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      # Check that the number of observation windows is consistent (upper coordinates)
      if(numObsWindows > dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      if(dim(lowerCoords)[2] != dimCoords | dim(upperCoords)[2] != dimCoords) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Perform the window reduction according to the local evaluation criterion
      redMatrix <- obsWindowReduce(lowerCoords, upperCoords, sourceCoords, localEvalParam, areAreas)
      numObsWindows <- dim(redMatrix)[1]
      if(numObsWindows <= 0) {
         # If there are no observation windows after the reduction process then return a zero likelihood
         if(log) {
            # Return the log scale density if requested
            return(-Inf)
         } else {
            # Return the natural scale density if requested
            return(0.0)
         }
      }
      inLowerCoords <- matrix(nrow = numObsWindows, ncol = dimCoords)
      inUpperCoords <- matrix(nrow = numObsWindows, ncol = dimCoords)
      inLowerCoords[1:numObsWindows, 1:dimCoords] <- redMatrix[1:numObsWindows, 1:dimCoords]                              # Retrieve the lower coordinates of the reduced windows
      inUpperCoords[1:numObsWindows, 1:dimCoords] <- redMatrix[1:numObsWindows, (dimCoords + 1):(dimCoords + dimCoords)]  # Retrieve the upper coordinates of the reduced windows
      # Recycle the intensity weights to match the number of observation windows
      recIntensityWeights <- numeric(length = numObsWindows, value = intensityWeights, recycle = TRUE)
      # Ensure that the intensity weights are valid
      if(sum(recIntensityWeights < 0.0) > 0) {
         # Invalid values for the intensity weights: set likelihood to zero
         if(log) {
            return(-Inf)
         } else {
            return(0.0)
         }
      }
      if(sum(recIntensityWeights) <= 0.0) {
         # If there are no points then this outcome is certain when the sum of the intensity
         # surface is zero
         if(numPoints == 0) {
            if(allowZero == 0.0) {
               if(log) {
                  return(-Inf)
               } else {
                  return(0.0)
               }
            } else {
               if(log) {
                  return(0.0)
               } else {
                  return(1.0)
               }
            }
            # If there is at least one point then the outcome is impossible when the sum of the
            # intensity surface is zero
         } else {
            if(log) {
               return(-Inf)
            } else {
               return(0.0)
            }
         }
      }
      # Test the validity of the source coordinates
      if(length(sourceCoords) != dimCoords) {
         stop("dimensioality of source coordinates are not consistent")
      }
      # Test the validity of the decay parameter
      if(normSD <= 0.0) {
         if(log) {
            # Return the log scale density if requested
            return(-Inf)
         } else {
            # Return the natural scale density if requested
            return(0.0)
         }
      }
      ## 5.1.3. Calculate the likelihood of each point ----
      # Calculate the integration of the decay kernel (weighted by intensity) for each observation window
      obsWindowNorm <- calcMNormSourceInt(
         matrix(inLowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         matrix(inUpperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         sourceCoords, normSD, areAreas) * recIntensityWeights
      totalNorm <- sum(obsWindowNorm)
      logPointDens <- rep(0.0, numPoints)
      # Initialise the output probability to the void probability
      outProb <- -totalNorm
      if(totalNorm <= 0.0) {
         # If there are no points then this outcome is certain when the sum of the intensity
         # surface is zero
         if(numPoints == 0) {
            if(allowZero == 0.0) {
               if(log) {
                  return(-Inf)
               } else {
                  return(0.0)
               }
            } else {
               if(log) {
                  return(0.0)
               } else {
                  return(1.0)
               }
            }
         # If there is at least one point then the outcome is impossible when the sum of the
         # intensity surface is zero
         } else {
            if(log) {
               return(-Inf)
            } else {
               return(0.0)
            }
         }
      } else if(numPoints > 0) {
         # Calculate the likelihood for each point
         for(pointIter in 1:numPoints) {
            curCoords <- x[pointIter, 1:dimCoords]
            # Retrieve the observation windows which the point falls within
            isInObsWindow <- isWithinWindow(curCoords,
               matrix(inLowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
               matrix(inUpperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
               tolDist, areAreas)
            # Calculate the distance decay component of the current point from the source (the exponent
            # of the isotropic multivariate normal diistribution)
            distDecay <- 0.0
            for(dimIter in 1:dimCoords) {
               distDecay <- distDecay + (curCoords[dimIter] - sourceCoords[dimIter]) * (curCoords[dimIter] - sourceCoords[dimIter])
            }
            distDecay <- distDecay * (-0.5 / (normSD * normSD))
            # Calculate the sum of the intensity for each point
            pointSumIntensity <- sum(isInObsWindow * recIntensityWeights * exp(distDecay))
            if(pointSumIntensity <= 0.0) {
               if(log) {
                  # Return the log scale density if requested
                  return(-Inf)
               } else {
                  # Return the natural scale density if requested
                  return(0.0)
               }
            }
            logPointDens[pointIter] <- log(pointSumIntensity)
         }
         outProb <- sum(logPointDens) - totalNorm# - sum(log(1:numPoints))
      } else if(allowZero == 0.0) {
         # Zero points are not possible when the distribution is truncated
         if(log) {
            return(-Inf)
         } else {
            return(0.0)
         }
      }
      ## 5.1.4. Correct the density for zero truncation (if appropriate) ----
      if(allowZero == 0.0) {
         outProb <- outProb - log(1.0 - exp(-totalNorm))
      }
      ## 5.1.5. Return the retrieved density ----
      if(log == 0) {
         # Export the output probability on the natural scale if requested
         outProb <- exp(outProb)
      }
      return(outProb)
   }
)

### 5.2. ==== Define the sampling function ====
# Define a function to draw random coordinates from the multivariate-normal source point inhomogenous Poisson process
rpoisMNormSourcePP <- nimbleFunction(
   run = function(
      n = integer(0),                              # Number of samples to draw from the distribution
      lowerCoords = double(2),                     # The lower coordinate values of the observation windows
      upperCoords = double(2),                     # The upper coordinate values of the observation windows
      sourceCoords = double(1),                    # The coordinates of the source location (the origin of the decay kernel)
      normSD = double(0),                          # The standard deviation of the isotropic multivariate normal distribution decay kernel
      intensityWeights = double(1, default = 1),   # Extra intensity weights for the different observation windows (by default a pure multivariate source point is assumed)
      areAreas = double(0, default = 1),           # Flag denoting whether the lower and upper coordinates are areas or transects
      numSamples = double(0, default = -1),        # Number of points (if negative then the number of rows in x is used to define this value)
      numWindows = double(0, default = -1),        # Number of observation windows (if negative the number of rows in lowerCoords is used to define this value)
      allowZero = double(0, default = 1),          # Allow zero samples? Non-zero values denote that there can be
      localEvalParam = double(0, default = -1)     # Parameter that controls the maximum distance that an observation window can be from x to be considered as a possible destination (-1 is interpredted as +Inf)
   ) {
      ## 5.2.1. Specify the return type dimensionality ----
      # Return type declaration
      returnType(double(2))
      ## 5.2.2. Sanity tests the inputs ----
      # Ensure that only one sample is requested
      if(n <= 0) { 
         stop("the number of requested samples must be above zero")
      } else if(n > 1) {
         print("rpoisMNormSourcePP only allows n = 1; using n = 1")
      }
      # Assess the dimensionality of the input coordinates
      dimCoords <- dim(lowerCoords)[2]
      # Ensure that the dimensionality is valid
      if(dimCoords <= 0) {
         stop("invalid dimension structure for the input coordinates")
      }
      # Assess the number of observation windows
      numObsWindows <- trunc(numWindows)
      if(numObsWindows <= 0) {
         numObsWindows <- dim(lowerCoords)[1]
      }
      # Ensure that the number of observation windows is valid
      if(numObsWindows <= 0) {
         stop("invalid number of observation windows")
      }
      # Check that the number of observation windows is consistent (lower coordinates)
      if(numObsWindows > dim(lowerCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      # Check that the number of observation windows is consistent (upper coordinates)
      if(numObsWindows > dim(upperCoords)[1]) {
         stop("number of observation windows not consistent between lower and upper coordinates (or the 'number of windows' parameter)")
      }
      if(dim(upperCoords)[2] != dimCoords) {
         stop("lower and/or upper coordinates have an incorrect dimension structure")
      }
      # Recycle the intensity weights to match the number of observation windows
      recIntensityWeights <- numeric(length = numObsWindows, value = intensityWeights, recycle = TRUE)
      # Ensure that the intensity weights are valid
      if(sum(recIntensityWeights < 0.0) > 0) {
         recIntensityWeights <- recIntensityWeights * numeric(value = recIntensityWeights > 0.0, length = numObsWindows)
      }
      # Weight the areas by their intensity weights multiplied by the integral of the decay function
      areaIntensityWeights <- recIntensityWeights * calcMNormSourceInt(
         matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         sourceCoords, normSD, areAreas)
      sumAreaIntensityWeights <- sum(areaIntensityWeights)
      if(sumAreaIntensityWeights <= 0.0) {
         return(matrix(nrow = 0, ncol = dimCoords, type = "double"))
      }
      ## 5.2.3. Generate a random number of points ----
      # Over-ride the random number of points if the number of samples has been set
      numPoints <- trunc(numSamples)
      if(numSamples < 0.0) {
         # If the number of samples has not been set then draw it from a Poisson distribution
         numPoints <- rpois(1, sumAreaIntensityWeights)
         if(allowZero == 0.0) {
            # If zeros are not possible then keep redrawing samples until a non-zero sample is drawn
            while(numPoints == 0) {
               numPoints <- rpois(1, sumAreaIntensityWeights)
            }
         }
      }
      if(numPoints < 0) {
         stop("invalid number of points to simulate")
      } else if(numPoints == 0) {
         if(allowZero == 0.0) {
            # Incompatible parameterisation
            stop("zero number of points requested for truncated distribution that excludes zero")
         }
         # Return an empty matrix (0 rows and dimCoords columns)
         return(matrix(nrow = 0, ncol = dimCoords, type = "double"))
      }
      ## 5.2.4. Generate random location for the points ----
      # Sample the points according to the multivariate-normal source point inhomogenous binomial process
      outCoordinates <- matrix(value = 0, nrow = numPoints, ncol = dimCoords)
      outCoordinates[1:numPoints, 1:dimCoords] <- rbinomMNormSourcePP(
         1, numPoints,
         matrix(lowerCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         matrix(upperCoords[1:numObsWindows, 1:dimCoords], nrow = numObsWindows, ncol = dimCoords),
         sourceCoords, normSD, intensityWeights, areAreas, numObsWindows, localEvalParam)
      return(outCoordinates)
   }
)

## 6. ------ REGISTER THE DISTRIBUTIONS ------
# Register the distributions with NIMBLE
registerDistributions(list(
   ### 6.1. ==== Register the dbinomPP distribution ====
   dbinomPP = list(
      ## 6.1.1. Define the BUGS code to call the distribution ----
      BUGSdist = "dbinomPP(numPoints, lowerCoords, upperCoords, intensityWeights, areAreas, numWindows)",
      ## 6.1.2. Set the input and output types and dimension structure ----
      types = c(
         "value = double(2)", "numPoints = double(0)", "lowerCoords = double(2)",
         "upperCoords = double(2)", "intensityWeights = double(1)", "areAreas = double(0)", "numWindows = double(0)"),
      ## 6.1.3. Define the cumulative probability and quantile function availability ----
      pqAvail = FALSE
   ),
   ### 6.2. ==== Register the dbinomPPSingle distribution ====
   dbinomPPSingle = list(
      ## 6.2.1. Define the BUGS code to call the distribution ----
      BUGSdist = "dbinomPPSingle(lowerCoords, upperCoords, intensityWeights, areAreas, numWindows)",
      ## 6.2.2. Set the input and output types and dimension structure ----
      types = c(
         "value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
         "intensityWeights = double(1)", "areAreas = double(0)", "numWindows = double(0)"),
      ## 6.2.3. Define the cumulative probability and quantile function availability ----
      pqAvail = FALSE
   ),
   ### 6.3. ==== Register the dbinomMNormSourcePP distribution ====
   dbinomMNormSourcePP = list(
      ## 6.3.1. Define the BUGS code to call the distribution ----
      BUGSdist = "dbinomMNormSourcePP(numPoints, lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, localEvalParam)",
      ## 6.3.2. Set the input and output types and dimension structure ----
      types = c(
         "value = double(2)", "numPoints = double(0)", "lowerCoords = double(2)",
         "upperCoords = double(2)", "sourceCoords = double(1)", "normSD = double(0)",
         "intensityWeights = double(1)", "areAreas = double(0)", "numWindows = double(0)",
         "localEvalParam = double(0)"),
      ## 6.3.3. Define the cumulative probability and quantile function availability ----
      pqAvail = FALSE
   ),
   ### 6.4. ==== Register the dbinomMNormSourcePPSingle distribution ====
   dbinomMNormSourcePPSingle = list(
      ## 6.4.1. Define the BUGS code to call the distribution ----
      BUGSdist = "dbinomMNormSourcePPSingle(lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, localEvalParam)",
      ## 6.4.2. Set the input and output types and dimension structure ----
      types = c(
         "value = double(1)", "lowerCoords = double(2)", "upperCoords = double(2)",
         "sourceCoords = double(1)", "normSD = double(0)", "intensityWeights = double(1)",
         "areAreas = double(0)", "numWindows = double(0)", "localEvalParam = double(0)"),
      ## 6.4.3. Define the cumulative probability and quantile function availability ----
      pqAvail = FALSE
   ),
   ### 6.5. ==== Register the dbinomMNormSourcePPMulti distribution ====
   dbinomMNormSourcePPMulti = list(
      ## 6.5.1. Define the BUGS code to call the distribution ----
      BUGSdist = "dbinomMNormSourcePPMulti(lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, localEvalParam)",
      ## 6.5.2. Set the input and output types and dimension structure ----
      types = c(
         "value = double(2)", "lowerCoords = double(2)", "upperCoords = double(2)",
         "sourceCoords = double(2)", "normSD = double(0)", "intensityWeights = double(1)",
         "areAreas = double(0)", "numWindows = double(0)", "localEvalParam = double(0)"),
      ## 6.5.3. Define the cumulative probability and quantile function availability ----
      pqAvail = FALSE
   ),
   ### 6.6. ==== Register the dbinomMNormSourcePP_dbinomPP_SingleMixture distribution ====
   dbinomMNormSourcePP_dbinomPP_SingleMixture = list(
      ## 6.6.1. Define the BUGS code to call the distribution ----
      BUGSdist = "dbinomMNormSourcePP_dbinomPP_SingleMixture(mixtureParam, lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numWindows, localEvalParam)",
      ## 6.6.2. Set the input and output types and dimension structure ----
      types = c(
         "value = double(1)", "mixtureParam = double(0)", "lowerCoords = double(2)", "upperCoords = double(2)",
         "sourceCoords = double(1)", "normSD = double(0)", "intensityWeights = double(1)",
         "areAreas = double(0)", "numWindows = double(0)", "localEvalParam = double(0)"),
      ## 6.6.3. Define the cumulative probability and quantile function availability ----
      pqAvail = FALSE
   ),
   ### 6.7. ==== Register dpoisPP distribution ====
   dpoisPP = list(
      ## 6.7.1. Define the BUGS code to call the distribution ----
      BUGSdist = "dpoisPP(lowerCoords, upperCoords, intensityWeights, areAreas, numSamples, numWindows, allowZero)",
      ## 6.7.2. Set the input and output types and dimension structure ----
      types = c(
         "value = double(2)", "lowerCoords = double(2)", "upperCoords = double(2)",
         "intensityWeights = double(1)", "areAreas = double(0)", "numSamples = double(0)",
         "numWindows = double(0)", "allowZero = double(0)"),
      ## 6.7.3. Define the cumulative probability and quantile function availability ----
      pqAvail = FALSE
   ),
   ### 6.8. ==== Register dpoisMNormSourcePP ====
   dpoisMNormSourcePP = list(
      ## 6.8.1. Define the BUGS code to call the distribution ----
      BUGSdist = "dpoisMNormSourcePP(lowerCoords, upperCoords, sourceCoords, normSD, intensityWeights, areAreas, numSamples, numWindows, allowZero, localEvalParam)",
      ## 6.8.2. Set the input and output types and dimension structure ----
      types = c(
         "value = double(2)", "lowerCoords = double(2)", "upperCoords = double(2)", "sourceCoords = double(1)",
         "normSD = double(0)", "intensityWeights = double(1)", "areAreas = double(0)", "numSamples = double(0)",
         "numWindows = double(0)", "allowZero = double(0)", "localEvalParam = double(0)"),
      ## 6.8.3. Define the cumlative probability and quantile function availability ----
      pqAvail = FALSE
   )
))