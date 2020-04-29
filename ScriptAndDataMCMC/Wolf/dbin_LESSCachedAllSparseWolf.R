#' @title Function to create a NIMBLE custom distribution for faster SCR model runs.
#'
#' @description
#' \code{dbin_LESS} returns the likelihood of a given individual & spatial binary detection history y[i,1:n.detectors] 
#' 
#' @param x \code{Vector} of length n.detectors containing detection/non detections 
#' @param pZero \code{Numeric} variable denoting the baseline detection parameter of the half-normal detection function.
#' @param sigma \code{Numeric} variable denoting the scale parameter of the half-normal detection function.
#' @param sxy A \code{Vector} of length 2 with individual activity center x and y coordinates.
#' @param d2 A \code{Vector}  of dimensions n.detectors with the square distances from the activity centers to each detector obtained from the \code{calculateDistance}.
#' @param maxDist A \code{Numeric} with the maximum distance to the AC where detections are considered possible. This applies a local evaluation of the state space (Milleret et al. 2018 Ecology and Evolution)
#' @param indicator A \code{Logical}; if == 0 no detections possible; saves time. 
#' @param log A \code{integer} indicates whether the log-probability or the probabilty should be returned (default to normal)
#'
#' @examples
#' y[i,1:n.detectors,t] ~ dbern_LESS(p0, sigma, d2[i,1:n.detectors,t], maxDist,  z[i,t]==2)

#### 1.Density function ####
dbin_LESSCachedAllSparseWolf <- nimbleFunction(run = function( x = double(1)
                                                           , sxy = double(1)
                                                           , sigma = double(0)
                                                           , nbDetections = double(0)
                                                           , yDets = double(1)
                                                           , detector.xy = double(2)
                                                           , trials = double(1)
                                                           , detectorIndex = double(2)
                                                           , nDetectorsLESS = double(1)
                                                           , ResizeFactor = double(0, default = 1)
                                                           , maxNBDets = double(0)
                                                           , habitatID = double(2)                      
                                                           , maxDist = double(0, default = 0.0)
                                                           , indicator = double(0, default = 1.0)
                                                           , z = double(0)
                                                           , p0State = double(2)
                                                           , idResponse = double(0)
                                                           , detCountries = double(1)
                                                           , detRoads = double(1)
                                                           , detTracks = double(1)
                                                           , detSnow = double(1)
                                                           , betaResponse = double(0)  
                                                           , betaTracks = double(0)
                                                           , betaRoads = double(0)
                                                           , betaSnow = double(0)
                                                           , log = integer(0, default = 0)){
   # RETURN TYPE DECLARATION
   returnType(double(0))
   
   ## CHECK INPUT DIMENSIONS
   nDetectors <- length(trials)

   ## SHORTCUT IF INDIVIDUAL IS NOT AVAILABLE FOR DETECTION
   if(indicator == 0){
      if(nbDetections == 0){
         if(log == 0) return(1.0)
         else return(0.0)
      }else{
         if(log == 0) return(0.0)
         else return(-Inf)
      }
   }

   ## GET DETECTOR INDEX FROM THE HABITATID MATRIX
   sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
   index <- detectorIndex[sxyID,1:nDetectorsLESS[sxyID]]
   ## GET NECESSARY INFO
   n.detectors <- length(index)
   maxDist_squared <- maxDist*maxDist

   ## RECREATE Y
   y <- nimNumeric(length = nDetectors, value = 0, init = TRUE)
   if(nbDetections > 0){
      for(j in 1:nbDetections){
         y[yDets[j]] <- x[j]
         ## check if a detection is out of the "detection window"
         if(sum(yDets[j]==index)==0){
            if(log == 0) return(0.0)
            else return(-Inf)
         }
      }
   }


   ## CALCULATE THE LIKELIHOOD
   alpha <- -1.0 / (2.0 * sigma * sigma)


   logProb <- 0.0
   count <- 1
   index1 <- c(index, 0) # so the count is not an issue

   for(j in 1:nDetectors){
      if(index1[count] == j){ # IF J IS EQUAL TO THE RELEVANT DETECTOR
         d2 <- pow(detector.xy[j,1] - sxy[1], 2) + pow(detector.xy[j,2] - sxy[2], 2)
         pZero <- ilogit(logit(p0State[detCountries[j], z]) +
                               betaResponse * idResponse +
                               betaRoads * detRoads[j] +
                               betaSnow * detSnow[j] +
                               betaTracks * detTracks[j])

         p <- pZero * exp(alpha * d2)
         logProb <- logProb + dbinom(y[j], prob = p, size = trials[j], log = TRUE)
         count <- count + 1
      }#else{logProb <- logProb + dbinom(y[j], prob = 0, size = trials[j], log = TRUE)}
   }
   if(log)return(logProb)
   return(exp(logProb))
   
})

#### 2.Sampling function ####


# rbin_LESSCachedAll<- nimbleFunction(run = function( n = integer(0)
#                                                     , pZero = double(0)
#                                                     , sigma = double(0)
#                                                     , detector.xy = double(2)
#                                                     , trials = double(1)
#                                                     , detectorIndex = double(2)
#                                                     , nDetectorsLESS = double(1)
#                                                     , ResizeFactor = double(0, default = 1)
#                                                     , maxNBDets = double(0)
#                                                     , habitatID = double(2)
#                                                    , maxDist = double(0, default = 0.0)
#                                                    , indicator = double(0, default = 1.0)){
#   # Return type declaration
#   returnType(double(1))
#   
#   # Check input dimensions
#   n.detectors <- length(d2)
#   
#   ## Shortcut if individual is not available for detection
#   if(indicator == 0){return(rep(0.0, n.detectors))}
#   
#   
#   ## Simulate the detections using the calculated detection distances ----
#   alpha <- -1.0 / (2.0 * sigma * sigma)
#   # Initialise a detections output vector
#   detectOut <- rep(0, n.detectors)
#   for(j in 1:n.detectors){
#     # Calculate the detection probability (negative distances repesent zero detection probability)
#     if(d2[j] <= (maxDist*maxDist)){
#       p <- pZero * exp(alpha * d2[j])
#       # Draw from a Bernoulli distribution with the calculated probability
#       detectOut[j] <- rbinom(1, trials[j], p)
#     }#if
#   }#j
#   
#   # Output
#   return(detectOut)
# })

#### 3.Registration ####
registerDistributions(list(
   dbin_LESSCachedAllSparseWolf = list(
      BUGSdist = "dbin_LESSCachedAllSparseWolf(sxy, sigma, nbDetections, yDets, 
      detector.xy, trials, detectorIndex, nDetectorsLESS, ResizeFactor, maxNBDets, habitatID, maxDist, indicator,
      z, p0State, idResponse, detCountries, detRoads, detTracks, detSnow, betaResponse, betaTracks, betaRoads, betaSnow)",
      types = c( "value = double(1)","sxy = double(1)", "sigma = double(0)",
                 "nbDetections = double(0)", "yDets = double(1)", "detector.xy = double(2)","trials = double(1)","detectorIndex = double(2)" ,
                 "nDetectorsLESS = double(1)", "ResizeFactor = double(0)",
                 "maxNBDets = double(0)","habitatID= double(2)",
                 "maxDist = double(0)" ,"indicator = double(0)"
                 , "z = double(0)", "p0State = double(2)"
                 , "idResponse = double(0)", "detCountries = double(1)"
                 , "detRoads = double(1)", "detTracks = double(1)"
                 , "detSnow = double(1)", "betaResponse = double(0)"  
                 , "betaTracks = double(0)", "betaRoads = double(0)"
                 , "betaSnow = double(0)"),
      pqAvail = FALSE)))

