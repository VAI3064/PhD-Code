mynormalizeBetweenArrays <- function (object, method = NULL, targets = NULL, 
    ...) 
{
    if (is.null(method)) {
        if (is(object, "matrix")) {
            method = "quantile"
        }
        else if (is(object, "EListRaw")) {
            method = "quantile"
        }
        else {
            method = "Aquantile"
        }
    }
    choices <- c("none", "scale", "quantile", "Aquantile", "Gquantile", 
        "Rquantile", "Tquantile", "vsn", "cyclicloess")
    method <- match.arg(method, choices)
    if (method == "vsn") 
        stop("vsn method no longer supported. Please use normalizeVSN instead.")
    if (is(object, "matrix")) {
        if (!(method %in% c("none", "scale", "quantile", "cyclicloess"))) 
            stop("method not applicable to matrix objects")
        return(switch(method, none = object, scale = normalizeMedianValues(object), 
            quantile = normalizeQuantiles(object, ...), cyclicloess = normalizeCyclicLoess(object, 
                method = cyclic.method, ...)))
    }
    if (is(object, "EListRaw")) {
        if (method == "cyclicloess") 
            object$E <- Recall(log2(object$E), method = method, 
                ...)
        else object$E <- log2(Recall(object$E, method = method, 
            ...))
        object <- new("EList", unclass(object))
        return(object)
    }
    if (is(object, "RGList")) 
        object <- myMA.RG2(object)
    if (is.null(object$M) || is.null(object$A)) 
        stop("object doesn't appear to be RGList or MAList object")
    switch(method, scale = {
        object$M <- normalizeMedianAbsValues(object$M)
        object$A <- normalizeMedianAbsValues(object$A)
    }, quantile = {
        narrays <- NCOL(object$M)
        Z <- normalizeQuantiles(cbind(object$A - object$M/2, 
            object$A + object$M/2), ...)
        G <- Z[, 1:narrays]
        R <- Z[, narrays + (1:narrays)]
        object$M <- R - G
        object$A <- (R + G)/2
    }, Aquantile = {
        object$A <- normalizeQuantiles(object$A, ...)
    }, Gquantile = {
        G <- object$A - object$M/2
        E <- normalizeQuantiles(G, ...) - G
        object$A <- object$A + E
    }, Rquantile = {
        R <- object$A + object$M/2
        E <- normalizeQuantiles(R, ...) - R
        object$A <- object$A + E
    }, Tquantile = {
        narrays <- NCOL(object$M)
        if (NCOL(targets) > 2) targets <- targets[, c("Cy3", 
            "Cy5")]
        targets <- as.vector(targets)
        Z <- cbind(object$A - object$M/2, object$A + object$M/2)
        for (u in unique(targets)) {
            j <- targets == u
            Z[, j] <- normalizeQuantiles(Z[, j], ...)
        }
        G <- Z[, 1:narrays]
        R <- Z[, narrays + (1:narrays)]
        object$M <- R - G
        object$A <- (R + G)/2
    })
    object
}
