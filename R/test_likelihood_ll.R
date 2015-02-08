two.state.trans.prob = function(forward.rate, backward.rate, elapsed.time){
  total.rate = forward.rate + backward.rate               
  
  return((matrix(c(rep(backward.rate,2),rep(forward.rate,2)),2,2) +
            matrix(c(forward.rate, -backward.rate, -forward.rate, backward.rate),2,2)*
            exp(-total.rate*elapsed.time))/total.rate)    
}


inRTwoStatePartLikehoods = function(my.tree, my.data, forward.rate, backward.rate){
  
  ## reorder the edges in the "pruningwise" order
  my.tree = reorder(my.tree, order = "pr")
  
  if (!("phylo" %in% class(my.tree)))
    stop("Error: object \"my.tree\" is not of class \"phylo\"")
  
  if (is.null(my.tree$edge.length))
    stop("Error: tree \" my.tree\" must have branch lengths.")
  
  ## reorder data on tips to match the order of the my.tree phylo object
  if (!is.null(names(my.data))) {
    if(!any(is.na(match(names(my.data), my.tree$tip.label)))){
      my.data = my.data[my.tree$tip.label]
    }else{
      warning('the names of argument "my.data" and the names of the tip labels
did not match: the former were ignored in the analysis.')
    }
  }
  
  ## prepare transition probability matrices (this of course can and should be done in C++ as well)
  ## prob.array = array(0, dim=c(2,2,length(my.tree$edge.length)))            
  
  ## for (i in 1:length(my.tree$edge.length)){
  ##  prob.array[,,i] = two.state.trans.prob(forward.rate, backward.rate, my.tree$edge.length[i])          
  ## }            
  
  return(twoStateSufficientStatistics(my.tree$edge, my.data, my.tree$edge.length, forward.rate, backward.rate))
}
