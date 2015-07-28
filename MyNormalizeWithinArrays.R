mynormalizeWithinArrays <- function (object, layout = object$printer, method = "loess", 
    weights = object$weights, span = 0.3, iterations = 4, controlspots = NULL, 
    df = 5, robust = "M", bc.method = "normexp", normexp.method='mle',offset = 0) 
{
    if (!is(object, "MAList")) 
        object <- myMA.RG2(object, bc.method = bc.method, offset = offset)
    choices <- c("none", "median", "loess", "printtiploess", 
        "composite", "control", "robustspline")
    method <- match.arg(method, choices)
    if (method == "none") 
        return(object)
    if (is.vector(object$M)) 
        object$M <- as.matrix(object$M)
    nprobes <- nrow(object$M)
    narrays <- ncol(object$M)
    if (!is.null(weights)) 
        weights <- asMatrixWeights(weights, dim = c(nprobes, 
            narrays))
    if (method == "median") {
        if (is.null(weights)) 
            for (j in 1:narrays) object$M[, j] <- object$M[, 
                j] - median(object$M[, j], na.rm = TRUE)
        else for (j in 1:narrays) object$M[, j] <- object$M[, 
            j] - weighted.median(object$M[, j], weights[, j], 
            na.rm = TRUE)
        return(object)
    }
    if (is.vector(object$A)) 
        object$A <- as.matrix(object$A)
    if (nrow(object$A) != nprobes) 
        stop("row dimension of A doesn't match M")
    if (ncol(object$A) != narrays) 
        stop("col dimension of A doesn't match M")
    switch(method, loess = {
        for (j in 1:narrays) {
            y <- object$M[, j]
            x <- object$A[, j]
            w <- weights[, j]
            object$M[, j] <- loessFit(y, x, w, span = span, iterations = iterations)$residuals
        }
    }, printtiploess = {
        if (is.null(layout)) stop("Layout argument not specified")
        ngr <- layout$ngrid.r
        ngc <- layout$ngrid.c
        nspots <- layout$nspot.r * layout$nspot.c
        nprobes2 <- ngr * ngc * nspots
        if (nprobes2 != nprobes) stop("printer layout information does not match M row dimension")
        for (j in 1:narrays) {
            spots <- 1:nspots
            for (gridr in 1:ngr) for (gridc in 1:ngc) {
                y <- object$M[spots, j]
                x <- object$A[spots, j]
                w <- weights[spots, j]
                object$M[spots, j] <- loessFit(y, x, w, span = span, 
                  iterations = iterations)$residuals
                spots <- spots + nspots
            }
        }
    }, composite = {
        if (is.null(layout)) stop("Layout argument not specified")
        if (is.null(controlspots)) stop("controlspots argument not specified")
        ntips <- layout$ngrid.r * layout$ngrid.c
        nspots <- layout$nspot.r * layout$nspot.c
        for (j in 1:narrays) {
            y <- object$M[, j]
            x <- object$A[, j]
            w <- weights[, j]
            f <- is.finite(y) & is.finite(x)
            if (!is.null(w)) f <- f & is.finite(w)
            y[!f] <- NA
            fit <- loess(y ~ x, weights = w, span = span, subset = controlspots, 
                na.action = na.exclude, degree = 0, surface = "direct", 
                family = "symmetric", trace.hat = "approximate", 
                iterations = iterations)
            alpha <- global <- y
            global[f] <- predict(fit, newdata = x[f])
            alpha[f] <- (rank(x[f]) - 1)/sum(f)
            spots <- 1:nspots
            for (tip in 1:ntips) {
                y <- object$M[spots, j]
                x <- object$A[spots, j]
                w <- weights[spots, j]
                local <- loessFit(y, x, w, span = span, iterations = iterations)$fitted
                object$M[spots, j] <- object$M[spots, j] - alpha[spots] * 
                  global[spots] - (1 - alpha[spots]) * local
                spots <- spots + nspots
            }
        }
    }, control = {
        if (is.null(controlspots)) stop("controlspots argument not specified")
        for (j in 1:narrays) {
            y <- object$M[, j]
            x <- object$A[, j]
            w <- weights[, j]
            f <- is.finite(y) & is.finite(x)
            if (!is.null(w)) f <- f & is.finite(w)
            y[!f] <- NA
            fit <- loess(y ~ x, weights = w, span = span, subset = controlspots, 
                na.action = na.exclude, degree = 1, surface = "direct", 
                family = "symmetric", trace.hat = "approximate", 
                iterations = iterations)
            y[f] <- y[f] - predict(fit, newdata = x[f])
            object$M[, j] <- y
        }
    }, robustspline = {
        if (is.null(layout)) stop("Layout argument not specified")
        for (j in 1:narrays) object$M[, j] <- normalizeRobustSpline(object$M[, 
            j], object$A[, j], layout, df = df, method = robust)
    })
    object
}
