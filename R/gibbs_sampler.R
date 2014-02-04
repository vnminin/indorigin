runTwoStateGibbs = function(inputTrees, rootDist, traitData, initLambda01, initLambda10, priorAlpha01, priorBeta01, priorAlpha10, priorBeta10, mcmcSize, mcmcBurnin, mcmcSubsample){
  
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
  
  mcmcOut = twoStatePhyloGibbsSampler(treeEdges, dim(treeEdges), treeBranchLengths, rootDist, treeTraits, 
                                initLambda01, initLambda10, priorAlpha01, priorBeta01, priorAlpha10, 
                                priorBeta10, mcmcSize, mcmcBurnin, mcmcSubsample)

  colnames(mcmcOut) = c("iter", "treeIndex", "logPost", "lambda01", "lambda10", "n01", "n10", "t0", "t1")
  
  
  return(list(treeList=inputTrees, mcmcOutput=mcmcOut))
  #return(treeEdges)
}