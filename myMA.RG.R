myMA.RG <- function (object, targets, bc.method = "none", offset = 0) 
{
    if (is.null(object$R) || is.null(object$G)) 
        stop("Object doesn't contain R and G components")
    #object <- backgroundCorrect(object, method = "subtract", offset = offset)
	targets<- as.matrix(targets)
    R <- object$R
    G <- object$G

  #  R[R <= 0] <- NA
  #  G[G <= 0] <- NA

    R <- log(R, 2)
    G <- log(G, 2)

 #  object$R <- object$G <- object$Rb <- object$Gb <- object$other <- NULL

for ( i in 1:length(targets[,2]))
	{
	if (targets[i,2] == "NG")
  		{
   	 	object$M <- as.matrix(R - G)
   	 	object$A <- as.matrix((R + G)/2)
		} 
	else
		{
	 	object$M <- as.matrix(G - R)
   	 	object$A <- as.matrix((G + R)/2)
		}
	}
new<- new("MAlist", unclass(object))
return(new)
}
