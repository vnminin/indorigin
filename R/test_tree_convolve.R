treeConvolveTest = function(my.tree, my.data, rate.01, rate.10,
		root.dist.prob.0, n.max) {
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

	res.prior = .Call("treeConvolve", my.tree$edge, my.tree$edge.length,
			my.data, rate.01, rate.10, n.max, root.dist.prob.0,
			"prior", TRUE,
			PACKAGE = "indorigin")
	res.posterior = .Call("treeConvolve", my.tree$edge, my.tree$edge.length,
            my.data, rate.01, rate.10, n.max, root.dist.prob.0,
            "posterior", TRUE,
            PACKAGE = "indorigin")
    res = cbind(res.prior, res.posterior)
	colnames(res) = c("prior", "posterior")
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