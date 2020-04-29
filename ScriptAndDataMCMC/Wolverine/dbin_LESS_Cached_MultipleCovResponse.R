#' @title NIMBLE custom distribution for faster SCR model runs.
#'
#' @description
#' \code{dbin_LESS_Cached_OneCov} returns the likelihood of a given individual & spatial binary detection history y[i,1:n.detectors] 
#' 
#' @param x \code{Vector} of length n.detectors containing detection/non detections 
#'
#' @examples
#' y[i,1:n.detectors,t] ~ dbern_LESS(p0, sigma, d2[i,1:n.detectors,t], maxDist,  z[i,t]==2)

#### 1.Density function ####
dbin_LESS_Cached_MultipleCovResponse <- nimbleFunction(run = function( x = double(1),
                                                          sxy = double(1),
                                                          sigma = double(0),
                                                          nbDetections = double(0),
                                                          yDets = double(1),
                                                          detector.xy = double(2),
                                                          trials = double(1),
                                                          detectorIndex = double(2),
                                                          nDetectorsLESS = double(1),
                                                          ResizeFactor = double(0, default = 1),
                                                          maxNBDets = double(0),
                                                          habitatID = double(2),                    
                                                          indicator = double(0, default = 1.0),
                                                          p0State = double(1),
                                                          detCountries = double(1),
                                                          detCov = double(2),
                                                          betaCov = double(1),
                                                          BetaResponse = double(0),
                                                          detResponse = double(0),
                                                          log = integer(0, default = 0)){
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
   
   ## GET DETECTOR INDEX FROM THE HABITAT ID MATRIX
   sxyID <- habitatID[trunc(sxy[2]/ResizeFactor)+1, trunc(sxy[1]/ResizeFactor)+1]
   index <- detectorIndex[sxyID,1:nDetectorsLESS[sxyID]]
   
   ## GET NECESSARY INFO 
   n.detectors <- length(index)
   n.covs <- dim(detCov)[2]
   alpha <- -1.0 / (2.0 * sigma * sigma)
   
   ## RECREATE Y 
   y <- nimNumeric(length = nDetectors, value = 0, init = TRUE)
   if(nbDetections > 0){
      for(j in 1:nbDetections){
         y[yDets[j]] <- x[j] 
         ## check if a detection is out of the "detection window"
         if(sum(yDets[j] == index) == 0){
            if(log == 0) return(0.0)
            else return(-Inf)
         } 
      }
   }
   
   ## CALCULATE THE LIKELIHOOD
   logProb <- 0.0 
   count <- 1 
   index1 <- c(index, 0) # so the count is not an issue
   
   for(j in 1:nDetectors){
      if(index1[count] == j){ # IF J IS EQUAL TO THE RELEVANT DETECTOR 
         d2 <- pow(detector.xy[j,1] - sxy[1], 2) + pow(detector.xy[j,2] - sxy[2], 2)
         pZero <- ilogit(logit(p0State[detCountries[j]]) +
                            sum(betaCov[1:n.covs] * detCov[j,1:n.covs]) +
                            BetaResponse*detResponse)
         p <- pZero * exp(alpha * d2)
         logProb <- logProb + dbinom(y[j], prob = p, size = trials[j], log = TRUE)
         count <- count + 1
      }
   }
   
   if(log)return(logProb)
   return(exp(logProb))
})



#### 3.Registration ####
registerDistributions(list(
   dbin_LESS_Cached_MultipleCovResponse = list(
      BUGSdist = "dbin_LESS_Cached_MultipleCovResponse(sxy, sigma, nbDetections, yDets, 
      detector.xy, trials, detectorIndex, nDetectorsLESS, ResizeFactor, maxNBDets,
      habitatID, indicator, p0State, detCountries, detCov, betaCov, BetaResponse, detResponse)",
      types = c( "value = double(1)","sxy = double(1)", "sigma = double(0)",
                 "nbDetections = double(0)", "yDets = double(1)", "detector.xy = double(2)",
                 "trials = double(1)", "detectorIndex = double(2)" ,
                 "nDetectorsLESS = double(1)", "ResizeFactor = double(0)",
                 "maxNBDets = double(0)","habitatID = double(2)",
                 "indicator = double(0)","p0State = double(1)",
                 "detCountries = double(1)", "detCov = double(2)", "betaCov = double(1)",
                 "BetaResponse = double(0)", "detResponse = double(0)"),
      pqAvail = FALSE)))
