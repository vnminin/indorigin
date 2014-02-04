testIndOrigin = function(inputTrees, rootDist=c(0,1), traitData, initLambda01, initLambda10, priorAlpha01, priorBeta01, priorAlpha10, priorBeta10, mcmcSize, mcmcBurnin, mcmcSubsample, mcSize, testThreshold=0){
  
  
  ## Prepare trees and trait data  for Rcpp code
  cat("pre-processing trees and trait data", "\n")
  
  if (!("multiPhylo" %in% class(inputTrees)))
    stop("Error: object \"inputTrees\" is not of class \"multiPhylo\"")
  
  treeNum = length(inputTrees)
  
  ## Allocate memory for an array to hold tree edge matrices
  treeEdges = array(0, dim=c(dim(inputTrees[[1]]$edge),treeNum))
  
  ## Allocate memory for a matrix to hold branch lengths
  treeBranchLengths = matrix(0, dim(inputTrees[[1]]$edge)[1], treeNum)
  
  ## Allocate memory for a matrix to hold data vectors (reordered for each tree)
  treeTraits = matrix(0, length(inputTrees[[1]]$tip.label), treeNum)
  
  for (i in 1:treeNum){
    
    if (is.null(inputTrees[[i]]$edge.length))
      stop("Error: All input trees must have branch lengths")
    
    if (!is.rooted(inputTrees[[i]]))
      stop("Error: All input trees must be rooted")
    
    ## reorder the edges in the "pruningwise" order
    tempTree = reorder(inputTrees[[i]], order = "pr")
    
    treeEdges[,,i] = tempTree$edge
    treeBranchLengths[,i] = tempTree$edge.length
    treeTraits[,i] = traitData
    
    ## reorder data on tips to match the order of the my.tree phylo object
    if (!is.null(names(traitData))){
      if(!any(is.na(match(names(traitData), tempTree$tip.label)))){
        treeTraits[,i] = traitData[tempTree$tip.label]
      }else{
        warning('the names of argument "inputTrees" and the names of the tip labels
did not match: the former were ignored in the analysis.')
      }
    }        
  }
  
  ## Run Gibbs sampler that iterates between drawing from the full conditional of missing data
  ## and drawing from the full conditional of model parameters (rates 0->1 and 1->0)
  
  cat("running Gibbs sampler", "\n")
  
  mcmcOut = twoStatePhyloGibbsSampler(treeEdges, dim(treeEdges), treeBranchLengths, rootDist, treeTraits, initLambda01, initLambda10, priorAlpha01, priorBeta01, priorAlpha10,
                                priorBeta10, mcmcSize, mcmcBurnin, mcmcSubsample)

  colnames(mcmcOut) = c("iter", "treeIndex", "logPost", "lambda01", "lambda10", "n01", "n10", "t0", "t1")
  
  cat("Computing posterior probabilities", "\n")
  
  ## compute the posterior probability Pr(N_01 <= testThreshold | data) 
  postProbs = numeric(dim(mcmcOut)[1])
  
  for (i in 1:length(postProbs)){
    postProbs[i] = sum(exp(treeConvolveTest(inputTrees[[mcmcOut[i,"treeIndex"]]], treeTraits[,mcmcOut[i,"treeIndex"]], mcmcOut[i,"lambda01"], mcmcOut[i,"lambda10"], 0, testThreshold)[,"posterior"]))
  }
  
  postProbEst = mean(postProbs)
  
  cat("Computing prior probabilities", "\n")  
  
  ## Produce a sample from the priors of lambda01 and lambda10
  priorSampleLambda01 = rgamma(mcSize, shape=priorAlpha01, rate=priorBeta01)
  priorSampleLambda10 = rgamma(mcSize, shape=priorAlpha10, rate=priorBeta10)
  
  ## compute the prior probability Pr(N_01 <= testThreshold) 
  priorProbs = numeric(mcSize)
  
  for (i in 1:mcSize){
    sampledTreeIndex = sample(c(1:treeNum),1)
    priorProbs[i] = sum(exp(treeConvolveTest(inputTrees[[sampledTreeIndex]], treeTraits[,sampledTreeIndex], priorSampleLambda01[i], priorSampleLambda10[i], 0, testThreshold)[,"prior"]))
  }
  
  priorProbEst = mean(priorProbs)
  
  bfVector = numeric(5)
  names(bfVector) = c(paste("Pr(N01<=",as.character(testThreshold),")",sep=""), paste("Pr(N01<=",as.character(testThreshold),"|data)",sep=""), paste("BF for N01<=",as.character(testThreshold),sep=""), "log10(BF)", "2xlog_e(BF)")
  bfVector[1] = priorProbEst
  bfVector[2] = postProbEst
  bfVector[3] = (postProbEst/(1-postProbEst))/(priorProbEst/(1-priorProbEst))
  bfVector[4] = log10(postProbEst) - log10(1-postProbEst) - log10(priorProbEst) + log10(1-priorProbEst)
  bfVector[5] = 2*(log(postProbEst) - log(1-postProbEst) - log(priorProbEst) + log(1-priorProbEst))  
  
  returnList = list(treeList=inputTrees, mcmcOutput=mcmcOut, priorJumpProbs=priorProbs,
                    postJumpProbs=postProbs, probsAndBFs=bfVector)
  
  class(returnList) = "indorigin"
  
  return(returnList)
}

getBF = function(indOriginResults){
  if (!("indorigin" %in% class(indOriginResults)))
    stop("Error: object \"indOriginResults\" is not of class \"indorigin\"")
  
  return(indOriginResults$probsAndBFs[3:5])
}

getPriorProb = function(indOriginResults){
  if (!("indorigin" %in% class(indOriginResults)))
    stop("Error: object \"indOriginResults\" is not of class \"indorigin\"")
  
  returnVector = numeric(3)
  names(returnVector) = c("PriorProbEst", "PriorProbSd", "PriorProbConfInt")
  
  est = indOriginResults$probsAndBFs[1]
  names(est) = NULL
  mcSd = sd(indOriginResults$priorJumpProbs)/sqrt(length(indOriginResults$priorJumpProbs))
  confInt = c(est-1.96*mcSd,est+1.96*mcSd)
  
  return(list(PriorProbEst=est, PriorProbSd=mcSd, PriorProbConfInt=confInt))
}

getPostProb = function(indOriginResults){
  if (!("indorigin" %in% class(indOriginResults)))
    stop("Error: object \"indOriginResults\" is not of class \"indorigin\"")
  
  est = indOriginResults$probsAndBFs[2]
  names(est) = NULL
  mcmcPostJumpProbs = coda::as.mcmc(indOriginResults$postJumpProbs)
  
  mcSd = sd(indOriginResults$postJumpProbs)/sqrt(coda::effectiveSize(mcmcPostJumpProbs))
  names(mcSd) = NULL
  confInt = c(est-1.96*mcSd,est+1.96*mcSd)
  names(confInt) = c("2.5%", "97.5%")
  
  return(list(PostProbEst=est, PostProbSd=mcSd, PostProbConfInt=confInt))
}

plotPosterior = function(indOriginResults, trace=FALSE){
  if (!("indorigin" %in% class(indOriginResults)))
    stop("Error: object \"indOriginResults\" is not of class \"indorigin\"")
  
mcmcTable = coda::mcmc(indOriginResults$mcmcOutput[,c("lambda01", "lambda10", "n01", "n10")], thin=as.integer(dim(indOriginResults$mcmcOutput)[1]/5000))
plot(mcmcTable, trace=trace)

}

AvereragePriorProbabilities = function(inputTrees, rootDist=c(0,1), traitData, initLambda01, initLambda10, testThreshold=0){
  
  ## Prepare trees and trait data  for Rcpp code
  cat("pre-processing trees and trait data", "\n")
  
  if (!("multiPhylo" %in% class(inputTrees)))
    stop("Error: object \"inputTrees\" is not of class \"multiPhylo\"")
  
  treeNum = length(inputTrees)
  
  ## Allocate memory for an array to hold tree edge matrices
  treeEdges = array(0, dim=c(dim(inputTrees[[1]]$edge),treeNum))
  
  ## Allocate memory for a matrix to hold branch lengths
  treeBranchLengths = matrix(0, dim(inputTrees[[1]]$edge)[1], treeNum)
  
  ## Allocate memory for a matrix to hold data vectors (reordered for each tree)
  treeTraits = matrix(0, length(inputTrees[[1]]$tip.label), treeNum)
  
  for (i in 1:treeNum){
    
    if (is.null(inputTrees[[i]]$edge.length))
      stop("Error: All input trees must have branch lengths")
    
    if (!is.rooted(inputTrees[[i]]))
      stop("Error: All input trees must be rooted")
    
    ## reorder the edges in the "pruningwise" order
    tempTree = reorder(inputTrees[[i]], order = "pr")
    
    treeEdges[,,i] = tempTree$edge
    treeBranchLengths[,i] = tempTree$edge.length
    treeTraits[,i] = traitData
    
    ## reorder data on tips to match the order of the my.tree phylo object
    if (!is.null(names(traitData))){
      if(!any(is.na(match(names(traitData), tempTree$tip.label)))){
        treeTraits[,i] = traitData[tempTree$tip.label]
      }else{
        warning('the names of argument "inputTrees" and the names of the tip labels
                did not match: the former were ignored in the analysis.')
      }
      }        
    }
  
  ## Run Gibbs sampler that iterates between drawing from the full conditional of missing data
  ## and drawing from the full conditional of model parameters (rates 0->1 and 1->0)
  
  cat("running Gibbs sampler", "\n")
  
  priorResults = AverageTreesConvolve(treeEdges, dim(treeEdges), treeBranchLengths, treeTraits, initLambda01, initLambda10, testThreshold, rootDist[1], "prior", 1)
  
  
  return(priorResults)
}

