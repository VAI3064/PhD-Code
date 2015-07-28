myMA.RG2 <- function (object, bc.method = "none", offset = 0) 
  {
    if (is.null(object$R) || is.null(object$G)) 
        stop("Object doesn't contain R and G components")
   # object <- backgroundCorrect(object, method = "subtract", offset = offset)
    R <- object$R
    G <- object$G
    R[R <= 0] <- NA
    G[G <= 0] <- NA
    R <- log(R, 2)
    G <- log(G, 2)
    object$R <- object$G <- object$Rb <- object$Gb <- object$other <- NULL
    object$M <- as.matrix(R - G)
    object$A <- as.matrix((R + G)/2)
    new("MAList", unclass(object))
}
