#' Gibbs sampler for a phylogenetic model with a binary trait.
#'
#' Runs a Gibbs sampler for a two-state CMTC 
#' on a set of phylogenetic trees, with data 
#' at the tips of the trees.
#' 
#' @param inputTrees \code{multiPhylo} object containing a list of trees.
#' @param rootDist Numeric vector of size two with probabilities of the 
#' two states at the root of the tree.
#' @param traitData Vector of trait values, coded as 0 and 1.
#' @param initLambda01 Initial value of the 0 -> 1 rate.
#' @param initLambda10 Initial value of the 1 -> 0 rate.
#' @param priorAlpha01 Shape parameter of the Gamma prior on the 0 -> 1 rate.
#' @param priorBeta01 Rate parameter of the Gamma prior on the 0 -> 1 rate.
#' @param priorAlpha10 Shape parameter of the Gamma prior on the 1 -> 0 rate.
#' @param priorBeta10 Rate parameter of the Gamma prior on the 0 -> 1 rate.
#' @param mcmcSize Total number of MCMC iterations. 
#' @param mcmcBurnin Number of initial iterations to ignore.
#' @param mcmcSubsample Integer controlling how many MCMC iterations to save. 
#' For example, \code{mcmcSubsample = 10} saves every tenth iteration after ignoreing the first \code{mcmcBurnin} iterations.
#'
#' @return List with two elements: the input trees and a matrix of MCMC output.
#' The columns of MCMC output matrix are Iteration number, Index of a tree 
#' that was sampled at this iteration, Log posterior, 
#' 0 -> 1 rate, 1 -> 0 rate, Number of 0 -> 1 jumps, Number of 1 -> 0 jumps,
#' Time spent in state 0, Time spent in state 1
runTwoStateGibbs <- function(inputTrees, rootDist, traitData, initLambda01, initLambda10, priorAlpha01, priorBeta01, priorAlpha10, priorBeta10, mcmcSize, mcmcBurnin, mcmcSubsample){
  
  ## Prepare trees and trait data  for Rcpp code
  cat("pre-processing trees and trait data", "\n")
  
  if (!("multiPhylo" %in% class(inputTrees)))
    stop("Error: object \"inputTrees\" is not of class \"multiPhylo\"")
  
  treeNum <- length(inputTrees)
  
  ## Allocate memory for an array to hold tree edge matrices
  treeEdges <- array(0, dim=c(dim(inputTrees[[1]]$edge),treeNum))
  
  ## Allocate memory for a matrix to hold branch lengths
  treeBranchLengths <- matrix(0, dim(inputTrees[[1]]$edge)[1], treeNum)
  
  ## Allocate memory for a matrix to hold data vectors (reordered for each tree)
  treeTraits <- matrix(0, length(inputTrees[[1]]$tip.label), treeNum)
  
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
        treeTraits[,i] <- traitData[tempTree$tip.label]
      }else{
        warning('the names of argument "inputTrees" and the names of the tip labels
did not match: the former were ignored in the analysis.')
      }
    }        
  }
  
  ## Run Gibbs sampler that iterates between drawing from the full conditional of missing data
  ## and drawing from the full conditional of model parameters (rates 0->1 and 1->0)
  
  cat("running Gibbs sampler", "\n")
  
  mcmcOut <- twoStatePhyloGibbsSampler(treeEdges, dim(treeEdges), treeBranchLengths, rootDist, treeTraits, initLambda01, initLambda10, priorAlpha01, priorBeta01, priorAlpha10, priorBeta10, mcmcSize, mcmcBurnin, mcmcSubsample)

  colnames(mcmcOut) <- c("iter", "treeIndex", "logPost", "lambda01", "lambda10", "n01", "n10", "t0", "t1")
  
  
  return(list(treeList=inputTrees, mcmcOutput=mcmcOut))
  #return(treeEdges)
}