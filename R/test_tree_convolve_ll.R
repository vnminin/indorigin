treeConvolveTest = function(my.tree, my.data, rate.01, rate.10,
		root.dist.prob.0, n.max, dist.type) {
	if (!("phylo" %in% class(my.tree)))
		stop("Error: object \"my.tree\" is not of class \"phylo\"")
	
	if (is.null(my.tree$edge.length))
		stop("Error: tree \" my.tree\" must have branch lengths.")
	
	if (!is.rooted(my.tree))
	stop("Error: The input tree must be rooted")
  
	## reorder data on tips to match the order of the my.tree phylo object
	my.tree = reorder(my.tree, order = "postorder")
	if (!is.null(names(my.data))) {
		if (!any(is.na(match(names(my.data), my.tree$tip.label)))) {
			my.data = my.data[my.tree$tip.label]
        } else {
			warning('the names of argument "my.data" and the names of the tip labels
							did not match: the former were ignored in the analysis.')
		}
    }

	res = .Call("treeConvolve", my.tree$edge, my.tree$edge.length,
			my.data, rate.01, rate.10, n.max, root.dist.prob.0,
			dist.type, TRUE,
			PACKAGE = "indorigin")
	return(res)
}
 
  
treeConvolve = function(my.tree, my.data, rate.01, rate.10,
		root.dist.prob.0, n.max = NULL, dist.type = "posterior", gains = TRUE) {
	if (!("phylo" %in% class(my.tree)))
		stop("Error: object \"my.tree\" is not of class \"phylo\"")
    
	if (is.null(my.tree$edge.length))
		stop("Error: tree \" my.tree\" must have branch lengths.")
	
	if (!is.rooted(my.tree))
	stop("Error: The input tree must be rooted")
	
	if (any(is.na(my.data)))
		stop("Missing tip data not yet supported")
	
	if (!all(my.data %in% c(0, 1)))
		stop("all tip labels must equal zero or one")
    
    if (length(root.dist.prob.0) != 1)
        stop("root.dist.prob.0 must be a scalar")
    
    if (!(dist.type %in% c("prior", "posterior")))
        stop("\"dist.type\" must be either prior or posterior")
    
    if (is.null(n.max))
        stop("must specify maximum number of jumps over tree")

	## reorder data on tips to match the order of the my.tree phylo object
	my.tree = reorder(my.tree, order = "postorder")
	if (!is.null(names(my.data))) {
		if (!any(is.na(match(names(my.data), my.tree$tip.label)))) {
			my.data = my.data[my.tree$tip.label]
		} else {
			warning('the names of argument "my.data" and the names of the tip labels
							did not match: the former were ignored in the analysis.')
		}
	}

	res = .Call("treeConvolve", my.tree$edge, my.tree$edge.length,
			my.data, rate.01, rate.10, n.max, root.dist.prob.0, dist.type,
			as.integer(gains), PACKAGE = "indorigin")
	return(res)
}

manyTreesConvolveTest = function(selected.trees, my.trees, my.data, rate.01, rate.10, root.dist.prob.0, n.max, dist.type, gains = TRUE) {
  for (i in 1:length(my.trees)){
    if (!("phylo" %in% class(my.trees[[i]])))
      stop("Error: object \"my.tree\" is not of class \"phylo\"")
    
    if (is.null(my.trees[[i]]$edge.length))
      stop("Error: tree \" my.tree\" must have branch lengths.")
    
    if (!is.rooted(my.trees[[i]]))
      stop("Error: The input tree must be rooted")
    
    ## reorder data on tips to match the order of the my.tree phylo object
    my.trees[[i]] = reorder(my.trees[[i]], order = "postorder")
  }
  
  if (!is.null(names(my.data))) {
    if (!any(is.na(match(names(my.data), my.trees[[1]]$tip.label)))) {
      my.data = my.data[my.trees[[1]]$tip.label]
    } else {
      warning('the names of argument "my.data" and the names of the tip labels
              did not match: the former were ignored in the analysis.')
    }
    }
  
  edges <- unlist(sapply(my.trees, "[", i = "edge"))
  branch.lengths <- matrix(unlist(sapply(my.trees, "[", i = "edge.length")), ncol = length(my.trees))
  res = .Call("manyTreesConvolve", selected.trees, edges, c(nrow(my.trees[[1]]$edge), ncol(my.trees[[1]]$edge), length(my.trees)), branch.lengths, my.data, rate.01, rate.10, n.max, root.dist.prob.0, dist.type, as.integer(gains), PACKAGE = "indorigin")
  return(res)
}
