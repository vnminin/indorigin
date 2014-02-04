pcltTest <- function(rate.01, rate.10, t, n.max) {
    rout <- rep(0.0, n.max + 1)
    rpar <- c(t, rate.01, rate.10)
    res <- .Call("pclt_test", rpar, rout)
    return(as.vector(res))
}